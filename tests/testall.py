import os
import unittest

class ArtsRun:
    """A baby arts run class for handling test cases"""
    def __init__(self,control_file):
        """Initialised by the control file name"""
        self.control_file=control_file
    def run(self):
        """Run the control file"""
        w,r,e=os.popen3('../src/arts '+self.control_file)
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
            
class MonteCarloTest(unittest.TestCase):
    """Testing the ARTS-MC algorithm"""
    MCrun=ArtsRun('simpleMC.arts')
    def test1(self):
        """Simple Monte Carlo test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running simpleMC.arts: '+self.MCrun.error
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
        
class MonteCarloTest2(unittest.TestCase):
    """Testing the ARTS-MC algorithm with incoming lookup"""
    MCrun=ArtsRun('simpleMCWithIncomingLookup.arts')
    def test1(self):
        """Simple Monte Carlo test (with incoming lookup) should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running simpleMCWithIncomingLookup.arts: '+self.MCrun.error
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
        
class DoitTest(unittest.TestCase):
    """Testing the ARTS-DOIT algorithm"""
    DOITrun=ArtsRun('simpleDOIT.arts')
    def test1(self):
        """Simple DOIT test should run with no errors"""
        self.DOITrun.run()
        assert self.DOITrun.error=='','Error running simpleDOIT.arts: '+self.DOITrun.error
    def test2(self):
        """Total radiance should be close to 255.5 K"""
        I=self.DOITrun.get_val('y')[0]
        assert abs(I-255.) < 1, 'I (='+str(I)+'K) is too far away from 255.5 K'
    def test3(self):
        """Polarization difference should be close to 1.0 K"""
        Q=self.DOITrun.get_val('y')[1]
        assert abs(Q-1.) < 0.2, 'Q (='+str(Q)+'K) is too far away from 1. K'
        
if __name__=='__main__':
    unittest.main()


