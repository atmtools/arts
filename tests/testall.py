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
        assert abs(I-201.8) < 4, 'I (='+str(I)+'K) is too far away from 201.8 K'
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=self.MCrun.get_val('y')[1]
        assert abs(Q-7.6) < 1, 'Q (='+str(Q)+'K) is too far away from 7.6 K'
        


if __name__=='__main__':
    unittest.main()


