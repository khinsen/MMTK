# Run the complete test suite
#
# Written by Konrad Hinsen.
#

import unittest

import basic_tests
import universe_tests
import particleproperty_tests
import pickle_tests
import energy_tests
import trajectory_tests
import normal_mode_tests
import subspace_tests
import enm_tests
import ewald_tests
import pi_energy_trajectory_tests

def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTests(basic_tests.suite())
    test_suite.addTests(universe_tests.suite())
    test_suite.addTests(particleproperty_tests.suite())
    test_suite.addTests(pickle_tests.suite())
    test_suite.addTests(energy_tests.suite())
    test_suite.addTests(normal_mode_tests.suite())
    test_suite.addTests(subspace_tests.suite())
    test_suite.addTests(trajectory_tests.suite())
    test_suite.addTests(enm_tests.suite())
    test_suite.addTests(ewald_tests.suite())
    test_suite.addTests(pi_energy_trajectory_tests.suite())
    return test_suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())

