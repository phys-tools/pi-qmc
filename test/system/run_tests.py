import unittest
import glob
import os
import imp

def buildTestSuite():
    suite = unittest.TestSuite()
    for testcase in glob.glob('*/test_*.py'):
        module = imp.load_source(testcase[:-3], testcase)
        suite.addTest(unittest.TestLoader().loadTestsFromModule(module))
    return suite

if __name__ == "__main__":
    print "\npi-qmc integration test suite\n"
    print "Running system tests through python unittest suite.\n"
    suite = buildTestSuite()
    result = unittest.TextTestRunner(verbosity=20).run(suite)
