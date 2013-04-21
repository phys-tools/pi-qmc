import unittest
import subprocess
import os
import pitools
import math


class SHOTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.chdir("sho")
        out = file("pi.log", "w")
        process = subprocess.Popen("../../../bin/pi-qmc", 
            stdout=subprocess.PIPE, stdin=subprocess.PIPE, 
            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        process.wait()
        cls.h5file = pitools.openFile()

    @classmethod
    def tearDownClass(cls):
        cls.h5file.close()
        os.chdir("..")

    def test_energy(self):
        energy = self.h5file.getScalar("thermo_energy")
        e, de = energy.getAverage()
        expect = 0.75/math.tanh(0.25), 0.041
        self.assertAlmostEqual(e, expect[0], delta=0.08, msg=
            'wrong total energy, expected %f but got %f' % (expect[0], e))
        self.assertAlmostEqual(de, expect[1], delta=0.01, msg=
            'wrong error for energy, expected %f but got %f' % (expect[1], de))
