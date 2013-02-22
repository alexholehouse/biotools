import unittest
import sequenceTools
from aaGroups import *

class TestSequenceToolsFunctions(unittest.TestCase):

    # Build manager object for all tests here
    def setUp(self):
        pass

    def do_test(self):
        print "No tests to run"


    def test_calc_isMotifPresent(self):
        
        

        self.assertEqual(True, sequenceTools.calc_isMotifPresent("AFGHIKLLKPLKET", ["K","L","L","K"]))
        self.assertEqual(True, sequenceTools.calc_isMotifPresent("AFGHIKLLKPLKET", ["G","F","A"]))
        self.assertEqual(False, sequenceTools.calc_isMotifPresent("AFGHIKLLKPLKET", ["K","A","L","K"]))
        self.assertEqual(True, sequenceTools.calc_isMotifPresent("EDEDEDPEDEDDE", [NEG,HYDROPHOBIC,NEG]))
        self.assertEqual(True, sequenceTools.calc_isMotifPresent("GEIAQLWDF", [NEG,AROMATIC,NEG]))
        
        print ""
        print "ok = " + str(sequenceTools.calc_isMotifPresent("QEGFAEGFVRALAE", [AROMATIC,"*",NEG]))
        print "ok = " + str(sequenceTools.calc_isMotifPresent("DDVYNYLFD", [AROMATIC,"*",NEG]))
        

    
                
