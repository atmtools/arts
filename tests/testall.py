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
        w,r,e=os.popen3('cd '+self.subdir+'; ../../src/arts '+self.control_file)
        self.output=r.read()
        self.error=e.read()
    def get_val(self,name):
        """get a Numeric or vector value from standard output. Always returns a list
        .  This is where numpy would be nice"""
        str_list=self.output[self.output.index('*'+name+'*'):].splitlines()[1].split()
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
        """Total radiance should be close to 276.523 K"""
        I=self.CSrun.get_val('y')[0]
        assert abs(I-276.523) < 0.01, 'I (='+str(I)+'K) is too far away from 276.523 K'

        
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


