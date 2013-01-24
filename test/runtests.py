#!/usr/bin/python

import unittest
import __init__ as test



suite = unittest.TestLoader().loadTestsFromTestCase(test.klean_test.TestKleanFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

suite = unittest.TestLoader().loadTestsFromTestCase(test.sequenceTools_test.TestSequenceToolsFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

