import unittest
import subprocess
import os
import pitools
import math


class MadelungTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.chdir("madelung")
        out = file("pi.log", "w")
        process = subprocess.Popen("pi3D", stdout=subprocess.PIPE,
            stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        process.wait()
        cls.h5file = pitools.openFile()

    @classmethod
    def tearDownClass(cls):
        cls.h5file.close()
        os.chdir("..")

    def test_coulomb_energy(self):
        energy = self.h5file.getScalar("coulomb_energy")
        e, de = energy.getAverage()
        expect = -1.747564594633182190636
        self.assertAlmostEqual(e, expect, delta=1e-5, msg=
            'wrong coulomb energy, expected %f but got %f' % (expect, e))

    def test_energy(self):
        energy = self.h5file.getScalar("thermo_energy")
        e, de = energy.getAverage()
        expect = -1.747564594633182190636
        self.assertAlmostEqual(e, expect, delta=1e-4, msg=
            'wrong total energy, expected %f but got %f' % (expect, e))
