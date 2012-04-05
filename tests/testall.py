import os
import unittest
import artsXML
import tempfile
import numpy
import subprocess
from numpy import array
import string
import random

class ArtsRun:
    """A baby arts run class for handling test cases"""
    def __init__(self,subdir,control_file):
        """Initialised by the control file name"""
        self.subdir=subdir
        self.control_file=control_file
    def run(self):
        """Run the control file"""
        print self.subdir;

        artsbin='../../src/arts'
        try:
            if os.environ['CMAKE_GENERATOR'] == 'Xcode':
                artsbin='../../src/' + os.environ['CONFIGURATION'] + '/arts'
        except KeyError:
            pass

        w,r,e=os.popen3('cd ' + self.subdir + '; '
                + artsbin + ' ' + self.control_file)
        self.output=r.read()
        self.error=e.read()
    def get_val(self,name):
        """get a Numeric or vector value from standard output. Always returns a list
        .  This is where numpy would be nice"""
        startindex=self.output.index('*'+name+'*')
        endindex=self.output[startindex+1:].index('*'+name+'*')
        str_list=self.output[startindex:startindex+endindex].split()[1:]
        #convert to list of floats
        val_list=[]
        for s in str_list:
            val_list.append(float(s))
        return val_list

class Testatm_fields_compactAddSpecies(unittest.TestCase):
    """Test WSM atm_fields_compactAddSpecies
    """

    def setUp(self):
        # import PyARTS here so that only this part fails if it's missing
        self.cfile = tempfile.NamedTemporaryFile(suffix=".arts", delete=True)
        self.gf3file = tempfile.NamedTemporaryFile(delete=True)
        self.gf4file = tempfile.NamedTemporaryFile(delete=True)
        self.newgf4file = tempfile.NamedTemporaryFile(delete=True)

        from PyARTS.arts_types import GriddedField as gf
        self.GriddedField = gf
        from PyARTS import arts_file_components as afc
        self.afc = afc
        from PyARTS.general import quotify
        code = afc.ArtsCodeGenerator(generic=False)
        code.add(afc.WSM("GriddedField3Create", "gf3"))
        code.add(afc.ReadXML("gf3", self.gf3file.name))
        code.add(afc.ReadXML("atm_fields_compact", self.gf4file.name))
        code.add(afc.WSM("atm_fields_compactAddSpecies",
                         "atm_fields_compact", quotify("D"), "gf3"))
        code.add(afc.WSM("output_file_formatSetAscii"))
        code.add(afc.WriteXML("atm_fields_compact", self.newgf4file.name))
        code.toFile(self.cfile)

    def runcfile(self):
        """Runs the control-file
        """
        p = subprocess.Popen("arts "+self.cfile.name, shell=True,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, close_fds=True)

        stderr = p.stderr.read()
        self.assertTrue(stderr=="", "Arts run failed with message: %s" % stderr)

    def check_merge(self, gf4_old, gf3, p=None):
        """Checks if two GFs merged give correct result.

        Basically a wrapper around a minimal ARTS c-file running
        atm_fields_compactAddSpecies, although the c-file itself is
        constant, the fields change from test to test

        Third arg 'p' can be passed an index that will be log-ed.
        """

        gf = self.GriddedField
        gf3.save(self.gf3file.file)
        gf4_old.save(self.gf4file.file)
        self.runcfile()
        newgf4 = gf.load(self.newgf4file.name)
        name = "".join(random.choice(string.letters) for _ in xrange(10))
        # construct reference case
        #print gf4_old.data.shape, gf3.data[numpy.newaxis, ...].shape
        #print gf3.axes, 'to', gf4_old.axes[1:]
        if p is not None: # take log10
            gf3.axes[p] = numpy.log10(gf3.axes[p])
            gf4_old.axes[p+1] = numpy.log10(gf4_old.axes[p+1])
        ref = gf(gf4_old.axes[0].tolist()+[name],
                 gf4_old.axes[1],
                 gf4_old.axes[2],
                 gf4_old.axes[3],
                 numpy.concatenate((gf4_old.data,
                                    gf3.interpolate(*gf4_old.axes[1:]).data[numpy.newaxis, ...]),
                                   0))
        if p is not None: # revert log10-op
            gf3.axes[p] = gf3.axes[p]**10.
            gf4_old.axes[p+1] = gf4_old.axes[p+1]**10.

        #print ref.data - newgf4.data
        #print (ref.data - newgf4.data).max()
        mdiff = abs(ref.data-newgf4.data)
        if numpy.isnan(mdiff).any():
            print gf3.axes, "to", gf4_old.axes[1:]
        self.assertTrue((mdiff.max()<1e-7),
                        "Difference between old and new gf4." \
                        "To be interpolated:\n%s\n" \
                        "Expected: \n%s\n" \
                        "Got: \n%s\n" \
                        "Diff: \n%s\n" \
                        "Max diff: %e" % \
                        (gf3.data, ref.data, newgf4.data, mdiff, mdiff.max()))
    def tearDown(self):
        self.cfile.close() # also deletes the file
        self.gf3file.close()
        self.gf4file.close()
        self.newgf4file.close()

    def test_3D_equalgrid(self):
        r = numpy.arange
        gf3 = self.GriddedField(r(1.0, 5.0), r(5.0), r(6.0), r(4*5*6.0).reshape((4, 5, 6)))
        gf4 = self.GriddedField(["A", "B", "C"], r(1.0, 5.0), r(5.0), r(6.0),
                                r(3*4*5*6.0).reshape((3, 4, 5, 6)))
        self.check_merge(gf4, gf3, p=0)


    def test_3D_diffgrid_regular(self):
        r = numpy.arange
        ax1a, ax2a, ax3a = r(2.0, 5.1, 1.0), r(3.0, 7.1, 1.0), r(4.0, 9.1, 1.0)
        ax1b, ax2b, ax3b = r(2.0, 4.0, 0.5), r(3.0, 6.0, 0.5), r(4.0, 8.0, 0.5)
        gf3 = self.GriddedField(ax1a, ax2a, ax3a,
                                r(ax1a.size*ax2a.size*ax3a.size, dtype='float').reshape((ax1a.size, ax2a.size, ax3a.size)))
        gf4 = self.GriddedField(["A"], ax1b, ax2b, ax3b,
                                r(ax1b.size*ax2b.size*ax3b.size, dtype='float').reshape((1, ax1b.size, ax2b.size, ax3b.size)))
        self.check_merge(gf4, gf3, p=0)

    def test_3D_diffgrid_irregular(self):
        r = numpy.arange
        a = numpy.array
        ax1a, ax2a, ax3a = a([0.1, 1., 2.]), a([0., 1., 2.]), a([0., 1., 2.])
        ax1b, ax2b, ax3b = a([0.5, 1.5]), a([0.5, 1.5]), a([0.5, 1.5])
        gf3 = self.GriddedField(ax1a, ax2a, ax3a,
                                numpy.random.rand(ax1a.size, ax2a.size, ax3a.size).round(0))
        gf4 = self.GriddedField(["A"], ax1b, ax2b, ax3b,
                                numpy.random.rand(1, ax1b.size, ax2b.size, ax3b.size).round(0))
        self.check_merge(gf4, gf3, p=0)

    def test_1D_equalgrid(self):
        r = numpy.arange
        ax = r(0.01, 1e3, 17.0)
        gf3 = self.GriddedField(ax, [0], [0],
                                numpy.random.rand(ax.size, 1, 1))
        gf4 = self.GriddedField(["A"], ax, [0], [0],
                                numpy.random.rand(1, ax.size, 1, 1))
        self.check_merge(gf4, gf3, p=0)

if __name__=='__main__':
    unittest.main()


