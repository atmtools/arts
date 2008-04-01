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
    Amsurun=ArtsRun('AMSUB', 'TestAMSUB.arts')
    def test1(self):
        """AMSU-B test should run with no errors"""
        self.Amsurun.run()
        assert self.Amsurun.error=='','Error running TestAMSUB.arts: '+self.Amsurun.error
    def test2(self):
        """Total radiance should be close to the reference values"""
        Iref=array([[206.88824183225, 216.462102351192, 234.626811582318, 219.085569613317, 243.954571395148, 246.871777050885, 237.843232796618, 227.615816348413, 233.192780885998, 233.867243790764],
            [244.79749182399, 256.505996100133, 272.583939216468, 258.752084766613, 281.042039386401, 280.441528688159, 274.001499556087, 267.456107354327, 271.089356203675, 273.164307098658],
            [247.158533568627, 242.485387635278, 244.090239824494, 247.548132738101, 246.473351207385, 246.438063580444, 243.073686162242, 242.381165257308, 245.158903280948, 246.551021706208],
            [261.74288161284, 258.609204059642, 257.912410892722, 260.568141125811, 263.888416168005, 263.915825749411, 257.488376761352, 261.483546383459, 260.261596991756, 263.254660630527],
            [271.81923920704, 270.814443356476, 271.065533077755, 273.420342268019, 276.511784358578, 275.828414664139, 270.032139531133, 275.640344375078, 271.832587541114, 276.137791923373]])

        I = artsXML.load("AMSUB/AMSUB.ybatch.xml.generated")
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


