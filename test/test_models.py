import models
import unittest
import numpy as np


def read_octave_file(file_name: str):
    """Read a data file coming from Octave to validate our code."""
    data = []
    with open(file_name, "r") as f:
        # Skip the first 5 lines of the file.
        for i in range(5):
            f.readline()

        data = f.readline().strip()
        data = data.split(" ")

        # Do we have complex numbers?
        if data[0].startswith('('):
            dc = data.copy()
            data = []
            # Complex numbers are saved as tuples (re, im).
            for x in dc:
                re, im = x.lstrip('(').rstrip(')').split(',')
                data.append(float(re) + (1j * float(im)))
        else:
            data = [float(x) for x in data]

    return data


class AmplifierModelTest(unittest.TestCase):
    """Test if the amplifier model behaves as expected."""

    def test_amplifier_linear(self):
        input_data = read_octave_file("test/data/amp_input_linear.csv")
        output_data = read_octave_file("test/data/amp_output_linear.csv")

        pa = models.Amplifier()
        pa.mode = 'linear'
        pa.gain = 2.5
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_atan(self):
        input_data = read_octave_file("test/data/amp_input_atan.csv")
        output_data = read_octave_file("test/data/amp_output_atan.csv")

        pa = models.Amplifier()
        pa.mode = 'atan'
        pa.gain = 2.8
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_tanh(self):
        input_data = read_octave_file("test/data/amp_input_tanh.csv")
        output_data = read_octave_file("test/data/amp_output_tanh.csv")

        pa = models.Amplifier()
        pa.mode = 'tanh'
        pa.gain = 3.2
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_poly3(self):
        input_data = read_octave_file("test/data/amp_input_poly3.csv")
        output_data = read_octave_file("test/data/amp_output_poly3.csv")

        pa = models.Amplifier()
        pa.mode = 'poly3'
        pa.gain = 1.6
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_poly3_pm(self):
        input_data = read_octave_file("test/data/amp_input_poly3_pm.csv")
        output_data = read_octave_file("test/data/amp_output_poly3_pm.csv")

        pa = models.Amplifier()
        pa.mode = 'poly3_pm'
        pa.gain = 4.1
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_poly5(self):
        input_data = read_octave_file("test/data/amp_input_poly5.csv")
        output_data = read_octave_file("test/data/amp_output_poly5.csv")

        pa = models.Amplifier()
        pa.mode = 'poly5'
        pa.gain = 2.6
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_limiter(self):
        input_data = read_octave_file("test/data/amp_input_limiter.csv")
        output_data = read_octave_file("test/data/amp_output_limiter.csv")

        pa = models.Amplifier()
        pa.mode = 'limiter'
        pa.gain = 1.8
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_softlimiter(self):
        input_data = read_octave_file("test/data/amp_input_softlimiter.csv")
        output_data = read_octave_file("test/data/amp_output_softlimiter.csv")

        pa = models.Amplifier()
        pa.mode = 'softlimiter'
        pa.gain = 2.0
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))

    def test_amplifier_noise_var(self):
        """Test that the method for setting the noise variance is working properly."""
        pa = models.Amplifier()
        pa.set_noise_var(10, 2, 13)

        self.assertAlmostEqual(1.377329576e-19, pa.noise_var)

    def test_amplifier_amplitude_db(self):
        """Test that the method for setting the noise variance is working properly."""
        pa = models.Amplifier()
        pa.set_maximum_output_power(20)

        self.assertAlmostEqual(2.236067977, pa.max_output_amplitude)

    def test_amplifier_avg_power(self):
        """Test that the method for setting the noise variance is working properly."""
        pa = models.Amplifier()
        pa.set_average_power(20, 10)

        self.assertAlmostEqual(0.316227766, pa.gain)
