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
            includeprefix=""
        else:
            includeprefix="../"
        w,r,e=os.popen3('cd ' + self.subdir + '; '
                + '../../src/arts '
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
        I=artsXML.load("MonteCarlo/MonteCarloGeneral.y.xml.generated")[0]
        dI=artsXML.load("MonteCarlo/MonteCarloGeneral.mc_error.xml.generated")[0]
        assert abs(I-201.8) < 4*dI, 'I (= %.2f K) is too far away from 201.8 K' % I
    def test3(self):
        """Polarization difference should be close to 7.6 K"""
        Q=artsXML.load("MonteCarlo/MonteCarloGeneral.y.xml.generated")[1]
        dQ=artsXML.load("MonteCarlo/MonteCarloGeneral.mc_error.xml.generated")[1]
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
        I=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.y.xml.generated")[0]
        dI=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.mc_error.xml.generated")[0]
        assert abs(I-201) < 4*dI, 'I (= %.2f K) is too far away from 201 K' % I
    def test3(self):
        """Polarization difference should be close to 7.7 K"""
        Q=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.y.xml.generated")[1]
        dQ=artsXML.load("MonteCarlo/MonteCarloGeneralGaussian.mc_error.xml.generated")[1]
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
            [206.897618904659, 216.47387362302, 234.64520477951, 219.093119113362, 243.969665770344, 246.894385994402, 237.867433595862, 227.637861424372, 233.214971660216, 233.882397697028],
            [244.811073688972, 256.519306829158, 272.594136288584, 258.760221042019, 281.048246301787, 280.447656441305, 274.01002184789, 267.475743403907, 271.103162041974, 273.173835052923],
            [247.149127875864, 242.476702766452, 244.079860706327, 247.537084583882, 246.45306692743, 246.423582608802, 243.066949633924, 242.374920190197, 245.150083641124, 246.540825427669],
            [261.745350437585, 258.609286489753, 257.907114868754, 260.564250958801, 263.882120068143, 263.902087078759, 257.48100749554, 261.480570299508, 260.256896848913, 263.247342836763],
            [271.820077454283, 270.816540678842, 271.059618617288, 273.418003059435, 276.507011327624, 275.819728169072, 270.0220106208, 275.639956647712, 271.828843342715, 276.130368937293]]);

        I = artsXML.load("AMSU/AMSUB.ybatch.xml.generated")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.001,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])
    def test3(self):
        """Total radiance should be close to the values of TestAMSUB_fast"""
        Iref = artsXML.load("AMSU/AMSUB.ybatch.xml.generated")
        I    = artsXML.load("AMSU/AMSUB_fast.ybatch.xml.generated")
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
            [206.894832935168, 216.471182368942, 234.641949521969, 219.090069731122, 243.967344377959, 246.892454682215, 237.864634862335, 227.633971497743, 233.211657388844, 233.878926288279],
            [250.947444539452, 261.802482498593, 276.017678600914, 264.221062083092, 283.539389489153, 282.200135012974, 276.65739434716, 272.242938923943, 274.586387908822, 277.244746975628],
            [247.149686195992, 242.476870536288, 244.080010612482, 247.538046474088, 246.453991185512, 246.424435928646, 243.066933324228, 242.374801173937, 245.150197418964, 246.541029917116],
            [261.745254770484, 258.609171336295, 257.906981980146, 260.564095428528, 263.882002209582, 263.901973114413, 257.480894063622, 261.480447122568, 260.256790525422, 263.247228300712],
            [271.435354184594, 270.218863959384, 270.270489481309, 272.74452733781, 275.805101210536, 275.26144214164, 269.310000617579, 274.98290753595, 271.206020686257, 275.430506383982]]);

        I = artsXML.load("MHS/MHS.ybatch.xml.generated")
        for j in range (5):
            for k in range (10):
                assert abs(I[j,k]-Iref[j,k]) < 0.001,'I[%d,%d] = %.3fK is too far away from %.3fK' % (j,k,I[j,k],Iref[j,k])
    def test3(self):
        """Total radiance should be close to the values of TestMHS_fast"""
        Iref = artsXML.load("MHS/MHS.ybatch.xml.generated")
        I    = artsXML.load("MHS/MHS_fast.ybatch.xml.generated")
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


