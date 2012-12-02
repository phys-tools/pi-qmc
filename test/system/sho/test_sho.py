import unittest
import subprocess
import os
import pitools
import math


class SHOTestCase(unittest.TestCase):
    def setUp(self):
        os.chdir("sho")
        out = file("pi.log", "w")
        process = subprocess.Popen("pi3D", stdout=subprocess.PIPE,
            stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        process.wait()
        self.h5file = pitools.openFile()

    def tearDown(self):
        self.h5file.close()
        os.remove("pimc.dat")

    def test_energy(self):
        energy = self.h5file.getScalar("thermo_energy")
        e, de = energy.getAverage()
        expect = 0.75/math.tanh(0.25), 0.041
        self.assertAlmostEqual(e, expect[0], delta=0.01, msg=
            'wrong total energy, expected %f but got %f' % (expect[0], e))
        self.assertAlmostEqual(de, expect[1], delta=0.005, msg=
            'wrong error for energy, expected %f but got %f' % (expect[1], de))
