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
        assert abs(I-201.8) < 4*dI, 'I (='+str(I)+'K) is too far away from 201.8 K'
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=artsXML.load("MonteCarlo/MonteCarloGeneral.y.xml.generated")[1]
        dQ=artsXML.load("MonteCarlo/MonteCarloGeneral.mc_error.xml.generated")[1]
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
        I=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.y.xml.generated")[0]
        dI=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.mc_error.xml.generated")[0]
        assert abs(I-201) < 4*dI, 'I (='+str(I)+'K) is too far away from 201 K'
    def test3(self):
        """Polarization difference should be close to 7.7 K"""
        Q=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.y.xml.generated")[1]
        dQ=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.mc_error.xml.generated")[1]
        assert abs(Q-7.6) < 4*Q, 'Q (='+str(Q)+'K) is too far away from 7.6 K'

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
        assert abs(I-199.5) < 4*dI, 'I (='+str(I)+'K) is too far away from 199.5 K'
        
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
        assert abs( max(I)-113.2 ) < 0.1, 'I (='+str(max(I))+'K) is too far away from 113.2 K'

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
        I=artsXML.load("ClearSky/ClearSky.y1.xml.generated")[0]
        assert abs(I-249.68) < 0.01, 'I (='+str(I)+'K) is too far away from 249.68 K'
    def test3(self):
        """Difference between on-the-fly and lookup table should be below 0.01 K"""
        I1=artsXML.load("ClearSky/ClearSky.y1.xml.generated")[0]
        I2=artsXML.load("ClearSky/ClearSky.y2.xml.generated")[0]
        assert abs(I2-I1) < 0.02, 'Discrepancy (='+str(I2-I1)+'K) is too far away from 0 K'


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
            [283.922845658604, 285.433266218633, 292.005816296646, 290.814629720198, 294.637902054818, 291.896109522697, 289.884748329969, 291.37431966744, 289.838401593895, 295.09875484347],
            [281.79488297592, 282.567024112108, 287.279234477674, 287.301869950055, 290.66009951157, 287.350589389212, 284.987841197655, 288.16545188872, 285.365970230601, 290.766906684748],
            [247.158533568626, 242.485387635278, 244.090239824494, 247.5481327381, 246.473351207385, 246.438063580444, 243.073686162242, 242.381165257309, 245.158903280948, 246.551021706207],
            [261.742918749543, 258.609206351036, 257.912411012791, 260.568150892205, 263.888416168249, 263.915825751693, 257.488376783408, 261.483547654897, 260.261597049341, 263.254661016021],
            [272.221778675292, 270.904055405926, 271.075284996257, 273.550056008027, 276.5132021352, 275.829729379469, 270.03635743901, 275.671275961664, 271.840735738719, 276.15583394648]]);

        I = artsXML.load("AMSU/AMSUB.ybatch.xml.generated")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.01, 'I['+str(j)+','+str(k)+'] = '+str(I[j,k])+' K is too far away from '+str(Iref[j,k])+' K'

class TestAbs(unittest.TestCase):
    """Testing the ARTS Absorption module"""
    Absrun=ArtsRun('Abs', 'TestAbs.arts')
    def test1(self):
        """Simple Absorption test should run with no errors"""
        self.Absrun.run()
        assert self.Absrun.error=='','Error running TestAbs.arts: '+self.Absrun.error

if __name__=='__main__':
    unittest.main()


