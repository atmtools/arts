import os
import unittest
import artsXML
from numpy import array

class ArtsRun:
    """A baby arts run class for handling test cases"""
    def __init__(self,subdir,control_file):
        """Initialised by the control file name"""
        self.subdir=subdir
        self.control_file=control_file
    def run(self):
        """Run the control file"""
        print self.subdir;
        if os.environ['TOPSRCDIR'][0] == '/':
            includeprefix=''
        else:
            includeprefix='../'

        artsbin='../../src/arts'
        try:
            if os.environ['CMAKE_GENERATOR'] == 'Xcode':
                artsbin='../../src/' + os.environ['BUILD_STYLE'] + '/arts'
        except KeyError:
            pass

        w,r,e=os.popen3('cd ' + self.subdir + '; '
                + artsbin + ' '
                + '-I' + includeprefix + os.environ['TOPSRCDIR'] + '/includes '
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
        """Total radiance should be close to 112.18 K"""
        I=artsXML.load("ClearSky/TestClearSky.y1.xml")[0]
        assert abs(I-112.18) < 0.01, 'I (='+str(I)+'K) is too far away from 112.18 K'
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

if __name__=='__main__':
    unittest.main()


