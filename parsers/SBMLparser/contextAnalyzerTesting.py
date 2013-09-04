# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 13:31:27 2013

@author: proto
"""
import sys
sys.path.insert(0, '../utils/')
import consoleCommands as console

import unittest

class ParametrizedTestCase(unittest.TestCase):
    """ TestCase classes that want to be parametrized should
        inherit from this class.
    """
    def __init__(self, methodName='runTest', param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_klass, param=None):
        """ Create a suite containing all tests taken from the given
            subclass, passing them the parameter 'param'.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, param=param))
        return suite
        
    
class TestOne(ParametrizedTestCase):
    '''
    Test for file correctness
    '''
    def test_correctness(self):
        fileName = self.param
        score = console.correctness(fileName)
        self.assertTrue(score)


if __name__ == "__main__":      
    suite = unittest.TestSuite()
    #ran = range(1,410)
    suite.addTest(ParametrizedTestCase.parametrize(TestOne, param='complex/output19.bnglreduced.bngl'))
    unittest.TextTestRunner(verbosity=2).run(suite)

