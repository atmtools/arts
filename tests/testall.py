import os
import unittest

def dotest(name):
    exit_status=os.system('../src/arts -r 03 '+name+'.arts >/dev/null')
    assert exit_status==0,'Error occurred while running '+name+'.arts\n'


class MonteCarloTest(unittest.TestCase):
    """Testing the ARTS-MC algorithm"""
    def testRunsWithNoErrors(self):
        """Simple Monte Carlo test should run with no errors"""
        dotest('simpleMC')
        

class DOIT1DTest(unittest.TestCase):
    """Testing the DOIT-1D algorithm"""
    def testRunsWithNoErrors(self):
        """Simple 1D DOIT test should run with no errors"""
        dotest('simpleDOIT1D')

class DOIT3DTest(unittest.TestCase):
    """Testing the DOIT-3D algorithm"""
    def testRunsWithNoErrors(self):
        """Simple 3D DOIT test should run with no errors"""
        dotest('simpleDOIT3D')

if __name__=='__main__':
    unittest.main()


