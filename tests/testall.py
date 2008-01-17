import os
import unittest

class ArtsRun:
    """A baby arts run class for handling test cases"""
    def __init__(self,subdir,control_file):
        """Initialised by the control file name"""
        self.subdir=subdir
        self.control_file=control_file
    def run(self):
        """Run the control file"""
        print self.subdir;
        w,r,e=os.popen3('cd ' + self.subdir + '; '
                + '../../src/arts '
                + '-I' + "../" + os.environ['TOPSRCDIR'] + '/includes '
                + self.control_file)
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
            
class TestMonteCarloSimple(unittest.TestCase):
    """Testing the ARTS-MC algorithm"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloSimple.arts')
    def test1(self):
        """Simple Monte Carlo test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloSimple.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201.8 K"""
        I=self.MCrun.get_val('y')[0]
        dI=self.MCrun.get_val('mc_error')[0]
        assert abs(I-201.8) < 4*dI, 'I (='+str(I)+'K) is too far away from 201.8 K'
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=self.MCrun.get_val('y')[1]
        dQ=self.MCrun.get_val('mc_error')[1]
        assert abs(Q-7.6) < 4*Q, 'Q (='+str(Q)+'K) is too far away from 7.6 K'
        
class TestMonteCarloWithIncomingLookup(unittest.TestCase):
    """Testing the ARTS-MC algorithm with incoming lookup"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloWithIncomingLookup.arts')
    def test1(self):
        """Simple Monte Carlo test (with incoming lookup) should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloWithIncomingLookup.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201.8 K"""
        I=self.MCrun.get_val('y')[0]
        dI=self.MCrun.get_val('mc_error')[0]
        assert abs(I-201.8) < 4*dI, 'I (='+str(I)+'K) is too far away from 201.8 K'
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=self.MCrun.get_val('y')[1]
        dQ=self.MCrun.get_val('mc_error')[1]
        assert abs(Q-7.6) < 4*Q, 'Q (='+str(Q)+'K) is too far away from 7.6 K'
        
class TestMonteCarloGeneral(unittest.TestCase):
    """Testing the MCGeneral algorithm"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloGeneral.arts')
    def test1(self):
        """MCGeneral test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloSimple.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201.8 K"""
        I=self.MCrun.get_val('y')[0]
        dI=self.MCrun.get_val('mc_error')[0]
        assert abs(I-201.8) < 4*dI, 'I (='+str(I)+'K) is too far away from 201.8 K'
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=self.MCrun.get_val('y')[1]
        dQ=self.MCrun.get_val('mc_error')[1]
        assert abs(Q-7.6) < 4*Q, 'Q (='+str(Q)+'K) is too far away from 7.6 K'
        
class TestMonteCarloGeneralGaussian(unittest.TestCase):
    """Testing the MCGeneral algorithm with a Gaussian antenna response"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloGeneralGaussian.arts')
    def test1(self):
        """MCGeneral (Gaussian Antenna) test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloGeneralGaussian.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201 K"""
        I=self.MCrun.get_val('y')[0]
        dI=self.MCrun.get_val('mc_error')[0]
        assert abs(I-201) < 4*dI, 'I (='+str(I)+'K) is too far away from 201 K'
    def test3(self):
        """Polarization difference should be close to 7.7 K"""
        Q=self.MCrun.get_val('y')[1]
        dQ=self.MCrun.get_val('mc_error')[1]
        assert abs(Q-7.6) < 4*Q, 'Q (='+str(Q)+'K) is too far away from 7.6 K'
        
class TestDOIT(unittest.TestCase):
    """Testing the ARTS-DOIT algorithm"""
    DOITrun=ArtsRun('DOIT', 'TestDOIT.arts')
    def test1(self):
        """Simple DOIT test should run with no errors"""
        self.DOITrun.run()
        assert self.DOITrun.error=='','Error running TestDOIT.arts: '+self.DOITrun.error
    def test2(self):
        """Total radiance should be close to 204.5 K"""
        I=self.DOITrun.get_val('y')[0]
        assert abs(I-204.5) < 1., 'I (='+str(I)+'K) is too far away from 204.5 K'
    def test3(self):
        """Polarization difference should be close to 7.2 K"""
        Q=self.DOITrun.get_val('y')[1]
        assert abs(Q-7.2) < 1., 'Q (='+str(Q)+'K) is too far away from 7.2 K'
        

class TestClearSky(unittest.TestCase):
    """Testing clear sky calculations"""
    CSrun=ArtsRun('ClearSky', 'TestClearSky.arts')
    def test1(self):
        """Simple clear sky test should run with no errors"""
        self.CSrun.run()
        assert self.CSrun.error=='','Error running TestClearSky.arts: '+self.CSrun.error
    def test2(self):
        """Total radiance should be close to 249.68 K"""
        I=self.CSrun.get_val('vector_1')[0]
        assert abs(I-249.68) < 0.01, 'I (='+str(I)+'K) is too far away from 249.68 K'
    def test3(self):
        """Difference between on-the-fly and lookup table should be below 0.01 K"""
        I1=self.CSrun.get_val('vector_1')[0]
        I2=self.CSrun.get_val('y')[0]
        assert abs(I2-I1) < 0.02, 'Discrepancy (='+str(I2-I1)+'K) is too far away from 0 K'

        
class TestAMSUB(unittest.TestCase):
    """Testing AMSU-B calculations"""
    Amsurun=ArtsRun('AMSUB', 'TestAMSUB.arts')
    def test1(self):
        """AMSU-B test should run with no errors"""
        self.Amsurun.run()
        assert self.Amsurun.error=='','Error running TestAMSUB.arts: '+self.Amsurun.error
    def test2(self):
        """Total radiance should be close to the reference values"""
        Iref=[206.908, 216.487, 234.643, 219.107, 243.97, 246.883, 237.87, 227.632, 233.211, 233.885, 244.824, 256.532, 272.593, 258.775, 281.049, 280.445, 274.015, 267.47, 271.101, 273.177, 247.158, 242.484, 244.089, 247.547, 246.472, 246.437, 243.072, 242.379, 245.158, 246.549, 261.742, 258.607, 257.91, 260.566, 263.887, 263.914, 257.486, 261.48, 260.261, 263.253, 271.819, 270.812, 271.062, 273.417, 276.509, 275.827, 270.029, 275.638, 271.831, 276.135]

        for k in range (len(Iref)):
            I=self.Amsurun.get_val('ybatch')[k]
            assert abs(I-Iref[k]) < 0.01, 'I (='+str(I)+'K) is too far away from '+str(Iref[k])+' K'

class TestAbs(unittest.TestCase):
    """Testing the ARTS Absorption module"""
    Absrun=ArtsRun('Abs', 'TestAbs.arts')
    def test1(self):
        """Simple Absorption test should run with no errors"""
        self.Absrun.run()
        assert self.Absrun.error=='','Error running TestAbs.arts: '+self.Absrun.error
#     def test2(self):
#         """Total radiance should be close to 204.5 K"""
#         I=self.Absrun.get_val('y')[0]
#         assert abs(I-204.5) < 1., 'I (='+str(I)+'K) is too far away from 204.5 K'
#     def test3(self):
#         """Polarization difference should be close to 7.2 K"""
#         Q=self.Absrun.get_val('y')[1]
#         assert abs(Q-7.2) < 1., 'Q (='+str(Q)+'K) is too far away from 7.2 K'

if __name__=='__main__':
    unittest.main()


