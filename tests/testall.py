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

class TestMonteCarloDataPrepare(unittest.TestCase):
    """Preparing data for ARTS-MC tests"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloDataPrepare.arts')
    def test1(self):
        """MCDataPrepare.arts should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloDataPrepare.arts: '+self.MCrun.error
            
class TestMonteCarloGeneral(unittest.TestCase):
    """Testing the MCGeneral algorithm"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloGeneral.arts')
    def test1(self):
        """MCGeneral test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloSimple.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201.8 K"""
        I=artsXML.load("MonteCarlo/TestMonteCarloGeneral.y.xml")[0]
        dI=artsXML.load("MonteCarlo/TestMonteCarloGeneral.mc_error.xml")[0]
        assert abs(I-201.8) < 4*dI, 'I (= %.2f K) is too far away from 201.8 K' % I
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=artsXML.load("MonteCarlo/TestMonteCarloGeneral.y.xml")[1]
        dQ=artsXML.load("MonteCarlo/TestMonteCarloGeneral.mc_error.xml")[1]
        assert abs(Q-7.6) < 4*dQ, 'Q (= %.2f K) is too far away from 7.6 K' % Q
        
class TestMonteCarloGeneralGaussian(unittest.TestCase):
    """Testing the MCGeneral algorithm with a Gaussian antenna response"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloGeneralGaussian.arts')
    def test1(self):
        """MCGeneral (Gaussian Antenna) test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloGeneralGaussian.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201 K"""
        I=artsXML.load("MonteCarlo/TestMonteCarloGeneralGaussian.y.xml")[0]
        dI=artsXML.load("MonteCarlo/TestMonteCarloGeneralGaussian.mc_error.xml")[0]
        assert abs(I-201) < 4*dI, 'I (= %.2f K) is too far away from 201 K' % I
    def test3(self):
        """Polarization difference should be close to 7.7 K"""
        Q=artsXML.load("MonteCarlo/TestMonteCarloGeneralGaussian.y.xml")[1]
        dQ=artsXML.load("MonteCarlo/TestMonteCarloGeneralGaussian.mc_error.xml")[1]
        assert abs(Q-7.6) < 4*dQ, 'Q (= %.2f K) is too far away from 7.6 K' % Q

class TestRteCalcMC(unittest.TestCase):
    """Testing the WSM RteCalcMC"""
    MCrun=ArtsRun('MonteCarlo', 'TestRteCalcMC.arts')
    def test1(self):
        """RteCalcMC test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running RteCalcMC.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 199.5 K"""
        I=artsXML.load("MonteCarlo/TestRteCalcMC.y.xml")[0]
        dI=artsXML.load("MonteCarlo/TestRteCalcMC.y_error.xml")[0]
        assert abs(I-199.5) < 4*dI, 'I (=%.2f K) is too far away from 199.5 K' % I


class TestOdinSMR(unittest.TestCase):
    """Testing OdinSMR calculation"""
    ODINrun=ArtsRun('OdinSMR', 'TestOdinSMR.arts')
    def test1(self):
        """TestOdinSMR.arts should run with no errors"""
        self.ODINrun.run()
        assert self.ODINrun.error=='','Error running TestOdinSMR.arts: '+self.ODINrun.error
    def test2(self):
        """Max radiance should be close to 113.2 K"""
        I=artsXML.load("OdinSMR/TestOdinSMR.y.xml")
#        assert abs( max(I)-113.2 ) < 0.1, 'I (=%.2f K) is too far away from 113.2 K' % (max(I))


class TestDOIT(unittest.TestCase):
    """Testing the ARTS-DOIT algorithm"""
    DOITrun=ArtsRun('DOIT', 'TestDOIT.arts')
    def test1(self):
        """Simple DOIT test should run with no errors"""
        self.DOITrun.run()
        assert self.DOITrun.error=='','Error running TestDOIT.arts: '+self.DOITrun.error
    def test2(self):
        """Total radiance should be close to 204.5 K"""
        I=artsXML.load("DOIT/TestDOIT.y.xml")[0]
        assert abs(I-204.5) < 1., 'I (='+str(I)+'K) is too far away from 204.5 K'
    def test3(self):
        """Polarization difference should be close to 7.2 K"""
        Q=artsXML.load("DOIT/TestDOIT.y.xml")[1]
        assert abs(Q-7.2) < 1., 'Q (=%.2f K) is too far away from 7.2 K' % Q
        

class TestClearSky(unittest.TestCase):
    """Testing clear sky calculations"""
    CSrun=ArtsRun('ClearSky', 'TestClearSky.arts')
    def test1(self):
        """Simple clear sky test should run with no errors"""
        self.CSrun.run()
        assert self.CSrun.error=='','Error running TestClearSky.arts: '+self.CSrun.error
    def test2(self):
        """Total radiance should be close to 112.36 K"""
        I=artsXML.load("ClearSky/TestClearSky.y1.xml")[0]
        assert abs(I-112.36) < 0.01, 'I (='+str(I)+'K) is too far away from 112.36 K'
    def test3(self):
        """Difference between on-the-fly and lookup table should be below 0.01 K"""
        I1=artsXML.load("ClearSky/TestClearSky.y1.xml")[0]
        I2=artsXML.load("ClearSky/TestClearSky.y2.xml")[0]
        assert abs(I2-I1) < 0.01, 'Discrepancy (=%.3f K) is too large' % I2-I1

class TestGroundBased(unittest.TestCase):
    """Testing GroundBased calculation"""
    GBrun=ArtsRun('GroundBased', 'TestGbased.arts')
    def test1(self):
        """TestGbased.arts should run with no errors"""
        self.GBrun.run()
        assert self.GBrun.error=='','Error running TestGbased.arts: '+self.GBrun.error

class TestAMSUB(unittest.TestCase):
    """Testing AMSU-B calculations"""
    Amsurun=ArtsRun('AMSU', 'TestAMSUB.arts')
    def test1(self):
        """AMSU-B test should run with no errors"""
        self.Amsurun.run()
        assert self.Amsurun.error=='','Error running TestAMSUB.arts: '+self.Amsurun.error
    def test2(self):
        """Total radiance should be close to the reference values"""
        Iref=array([
            [206.9169, 216.4922, 234.6898, 219.1338, 244.0035, 246.9210, 237.9041, 227.6869, 233.2804, 233.9388],
            [244.8373, 256.5383, 272.6222, 258.8054, 281.0623, 280.4536, 274.0261, 267.5187, 271.1467, 273.2158],
            [247.1439, 242.4708, 244.0733, 247.5317, 246.4463, 246.4177, 243.0605, 242.3689, 245.1436, 246.5355],
            [261.7384, 258.6013, 257.8965, 260.5558, 263.8710, 263.8911, 257.4703, 261.4709, 260.2463, 263.2372],
            [271.8140, 270.8092, 271.0481, 273.4095, 276.4952, 275.8080, 270.0098, 275.6307, 271.8183, 276.1186]]);

        I = artsXML.load("AMSU/TestAMSUB.ybatch.xml")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.01,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])
    def test3(self):
        """Total radiance should be close to the values of TestAMSUB_fast"""
        Iref = artsXML.load("AMSU/TestAMSUB.ybatch.xml")
        I    = artsXML.load("AMSU/TestAMSUB_fast.ybatch.xml")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.2,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])

class TestAMSUB_fast(unittest.TestCase):
    """Testing AMSU-B fast calculation (reduced frequency set)"""
    Amsurun=ArtsRun('AMSU', 'TestAMSUB_fast.arts')
    def test1(self):
        """AMSU-B fast test should run with no errors"""
        self.Amsurun.run()
        assert self.Amsurun.error=='','Error running TestAMSUB_fast.arts: '+self.Amsurun.error

class TestMHS(unittest.TestCase):
    """Testing MHS calculations"""
    Mhsrun=ArtsRun('MHS', 'TestMHS.arts')
    def test1(self):
        """MHS test should run with no errors"""
        self.Mhsrun.run()
        assert self.Mhsrun.error=='','Error running TestMHS.arts: '+self.Mhsrun.error
    def test2(self):
        """Total radiance should be close to the reference values"""
        Iref=array([
            [206.9141, 216.4894, 234.6865, 219.1307, 244.0012, 246.9191, 237.9013, 227.6830, 233.2771, 233.9353],
            [250.9726, 261.8183, 276.0355, 264.2595, 283.5446, 282.1988, 276.6646, 272.2757, 274.6164, 277.2745],
            [247.1424, 242.4698, 244.0724, 247.5299, 246.4444, 246.4159, 243.0598, 242.3682, 245.1427, 246.5343],
            [261.7384, 258.6013, 257.8964, 260.5557, 263.8710, 263.8910, 257.4703, 261.4709, 260.2462, 263.2372],
            [271.4295, 270.2120, 270.2595, 272.7365, 275.7939, 275.2503, 269.2983, 274.9742, 271.1960, 275.4193]]);

        I = artsXML.load("MHS/TestMHS.ybatch.xml")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.01,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])
    def test3(self):
        """Total radiance should be close to the values of TestMHS_fast"""
        Iref = artsXML.load("MHS/TestMHS.ybatch.xml")
        I    = artsXML.load("MHS/TestMHS_fast.ybatch.xml")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.2,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])

class TestMHS_fast(unittest.TestCase):
    """Testing MHS fast calculation (reduced frequency set)"""
    Mhsrun=ArtsRun('MHS', 'TestMHS_fast.arts')
    def test1(self):
        """MHS fast test should run with no errors"""
        self.Mhsrun.run()
        assert self.Mhsrun.error=='','Error running TestMHS_fast.arts: '+self.Mhsrun.error

class TestAbs(unittest.TestCase):
    """Testing the ARTS Absorption module"""
    Absrun=ArtsRun('Abs', 'TestAbs.arts')
    def test1(self):
        """Simple Absorption test should run with no errors"""
        self.Absrun.run()
        assert self.Absrun.error=='','Error running TestAbs.arts: '+self.Absrun.error

class TestDOITBatch(unittest.TestCase):
    """Testing DOITBatch calculations"""
    Doitbatchrun=ArtsRun('DOITBatch', 'TestDOITBatch.arts')
    def test1(self):
        """DOITBatch test should run with no errors"""
        self.Doitbatchrun.run()
        assert self.Doitbatchrun.error=='','Error running TestDOITBatch.arts: '+self.Doitbatchrun.error
    def test2(self):
        """Total BT should be close to the reference values"""
        Iref=array([
                [244.3848,  243.2521,  270.6285,  229.2447,  282.4720,  282.6709],
                [176.9833,  176.9405,  271.1160,  259.5076,  261.2708,  262.5003],
                [259.0476,  257.7898,  258.9736,  218.6711,  275.5341,  273.9887],
                [204.4788,  204.0665,  277.7738,  257.0578,  268.1884,  269.2191]]);

        I = artsXML.load("DOITBatch/TestDOITBatch.ybatch.xml")
        for j in range (4):
            for k in range (6):
                assert abs(I[j,k]-Iref[j,k]) < 0.01,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])

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

class TestWfuns(unittest.TestCase):
    """Testing weighting function calculations"""
    TestWfunsrun=ArtsRun('Wfuns', 'TestWfuns.arts')
    def test1(self):
        """Wfuns test should run with no errors"""
        self.TestWfunsrun.run()
        assert self.TestWfunsrun.error=='','Error running TestWfuns.arts: '+self.TestWfunsrun.error

if __name__=='__main__':
    unittest.main()


