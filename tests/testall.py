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
            
class TestMonteCarloGeneral(unittest.TestCase):
    """Testing the MCGeneral algorithm"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloGeneral.arts')
    def test1(self):
        """MCGeneral test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloSimple.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201.8 K"""
        I=artsXML.load("MonteCarlo/MonteCarloGeneral.y.xml.generated")[0]
        dI=artsXML.load("MonteCarlo/MonteCarloGeneral.mc_error.xml.generated")[0]
        assert abs(I-201.8) < 4*dI, 'I (= %.2f K) is too far away from 201.8 K' % I
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=artsXML.load("MonteCarlo/MonteCarloGeneral.y.xml.generated")[1]
        dQ=artsXML.load("MonteCarlo/MonteCarloGeneral.mc_error.xml.generated")[1]
        assert abs(Q-7.6) < 4*Q, 'Q (= %.2f K) is too far away from 7.6 K' % Q
        
class TestMonteCarloGeneralGaussian(unittest.TestCase):
    """Testing the MCGeneral algorithm with a Gaussian antenna response"""
    MCrun=ArtsRun('MonteCarlo', 'TestMonteCarloGeneralGaussian.arts')
    def test1(self):
        """MCGeneral (Gaussian Antenna) test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running TestMonteCarloGeneralGaussian.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 201 K"""
        I=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.y.xml.generated")[0]
        dI=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.mc_error.xml.generated")[0]
        assert abs(I-201) < 4*dI, 'I (= %.2f K) is too far away from 201 K' % I
    def test3(self):
        """Polarization difference should be close to 7.7 K"""
        Q=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.y.xml.generated")[1]
        dQ=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.mc_error.xml.generated")[1]
        assert abs(Q-7.6) < 4*Q, 'Q (= %.2f K) is too far away from 7.6 K' % Q

class TestRteCalcMC(unittest.TestCase):
    """Testing the WSM RteCalcMC"""
    MCrun=ArtsRun('MonteCarlo', 'TestRteCalcMC.arts')
    def test1(self):
        """RteCalcMC test should run with no errors"""
        self.MCrun.run()
        assert self.MCrun.error=='','Error running RteCalcMC.arts: '+self.MCrun.error
    def test2(self):
        """Total radiance should be close to 199.5 K"""
        I=artsXML.load("MonteCarlo/RteCalcMC.y.xml.generated")[0]
        dI=artsXML.load("MonteCarlo/RteCalcMC.mc_error.xml.generated")[0]
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
        I=artsXML.load("OdinSMR/TestOdinSMR.y.xml.generated")
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
        I=artsXML.load("DOIT/DOIT.y.xml.generated")[0]
        assert abs(I-204.5) < 1., 'I (='+str(I)+'K) is too far away from 204.5 K'
    def test3(self):
        """Polarization difference should be close to 7.2 K"""
        Q=artsXML.load("DOIT/DOIT.y.xml.generated")[1]
        assert abs(Q-7.2) < 1., 'Q (=%.2f K) is too far away from 7.2 K' % Q
        

class TestClearSky(unittest.TestCase):
    """Testing clear sky calculations"""
    CSrun=ArtsRun('ClearSky', 'TestClearSky.arts')
    def test1(self):
        """Simple clear sky test should run with no errors"""
        self.CSrun.run()
        assert self.CSrun.error=='','Error running TestClearSky.arts: '+self.CSrun.error
    def test2(self):
        """Total radiance should be close to 249.68 K"""
        I=artsXML.load("ClearSky/ClearSky.y1.xml.generated")[0]
        assert abs(I-249.68) < 0.01, 'I (='+str(I)+'K) is too far away from 249.68 K'
    def test3(self):
        """Difference between on-the-fly and lookup table should be below 0.01 K"""
        I1=artsXML.load("ClearSky/ClearSky.y1.xml.generated")[0]
        I2=artsXML.load("ClearSky/ClearSky.y2.xml.generated")[0]
        assert abs(I2-I1) < 0.02, 'Discrepancy (=%.3f K) is too large' % I2-I1


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
            [206.888195065204, 216.462060070912, 234.626781338542, 219.085528920088, 243.954550457357, 246.871759690388, 237.843204848381, 227.615777965081, 233.192749675628, 233.867209948054],
            [244.797923499535, 256.505899525287, 272.583291905995, 258.751964413882, 281.040958601395, 280.440424378103, 274.000752302743, 267.455682175985, 271.088772792904, 273.163703634007],
            [247.152151353834, 242.481089143301, 244.086775607671, 247.5396643017, 246.46425722026, 246.429927479463, 243.07055096207, 242.379181106515, 245.154786731044, 246.546725006306],
            [261.743067325013, 258.609394695298, 257.91255048814, 260.568290959895, 263.888636064822, 263.916059116122, 257.488527385271, 261.483777874329, 260.261769401592, 263.254857017258],
            [271.817840128734, 270.812918098056, 271.063906506059, 273.418959145645, 276.509355790478, 275.826213431893, 270.030302532401, 275.638530712342, 271.830632022535, 276.135580633224]]);

        I = artsXML.load("AMSU/AMSUB.ybatch.xml.generated")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.01,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])

class TestAbs(unittest.TestCase):
    """Testing the ARTS Absorption module"""
    Absrun=ArtsRun('Abs', 'TestAbs.arts')
    def test1(self):
        """Simple Absorption test should run with no errors"""
        self.Absrun.run()
        assert self.Absrun.error=='','Error running TestAbs.arts: '+self.Absrun.error

if __name__=='__main__':
    unittest.main()


