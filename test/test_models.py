import models
import unittest
import numpy as np
import numpy.random
import models.utils


def read_octave_file(file_name: str):
    """Read a data file coming from Octave to validate our code."""
    data = []
    with open(file_name, "r") as f:
        # Skip the first 5 lines of the file.
        for i in range(5):
            f.readline()

        for line in f:
            line = line.strip()
            line = line.split(" ")

            # Do we have complex numbers?
            if line[0].startswith('('):
                ldata = []
                # Complex numbers are saved as tuples (re, im).
                for x in line:
                    re, im = x.lstrip('(').rstrip(')').split(',')
                    ldata.append(float(re) + (1j * float(im)))
                data.append(ldata)
            elif line[0] != '':
                ldata = [float(x) for x in line]
                data.append(ldata)

    return np.array(data).T


class CouplerModelTest(unittest.TestCase):
    """Test if the coupler model behaves as expected."""

    def test_coupler_0db(self):
        input_data = read_octave_file("test/data/coupler_input_10db.csv")

        coup = models.Coupler()
        coup.damping = 0
        coupler_data = coup.run(input_data)

        # With 0dB damping the input and output data should be the same.
        self.assertTrue(np.allclose(coupler_data, input_data))

    def test_coupler_10db(self):
        input_data = read_octave_file("test/data/coupler_input_10db.csv")
        output_data = read_octave_file("test/data/coupler_output_10db.csv")

        coup = models.Coupler()
        coup.damping = -10
        coupler_data = coup.run(input_data)

        self.assertTrue(np.allclose(coupler_data, output_data))


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

    def test_amplifier_noise_var(self):
        """Test if the noise variance works correctly."""
        numpy.random.seed(45612)
        input_data = read_octave_file("test/data/amp_input_noise.csv")
        output_data = read_octave_file("test/data/amp_output_noise.csv")

        pa = models.Amplifier()
        pa.mode = 'ideal'
        pa.set_noise_var(22, 2, 5)
        amp_data = pa.run(input_data)

        self.assertTrue(np.allclose(amp_data, output_data))


class DacModelTest(unittest.TestCase):
    """Test if the DAC model behaves as expected."""

    def test_dac(self):
        input_data = read_octave_file("test/data/dac_input.csv")
        output_data = read_octave_file("test/data/dac_output.csv")

        dac = models.Dac()
        dac.nobits = 10
        dac.trunc_level = 0.6
        dac_data = dac.run(input_data)

        self.assertTrue(np.allclose(dac_data, output_data))

    def test_dac_30nobits(self):
        input_data = read_octave_file("test/data/dac_input_30nobits.csv")
        output_data = read_octave_file("test/data/dac_output_30nobits.csv")

        dac = models.Dac()
        dac.nobits = 30
        dac.trunc_level = 0.58
        dac_data = dac.run(input_data)

        self.assertTrue(np.allclose(dac_data, output_data))

    def test_dac_inf(self):
        input_data = read_octave_file("test/data/dac_input_30nobits.csv")

        dac = models.Dac()
        dac_data = dac.run(input_data)

        # With an infinite limit the input and output data should be the same.
        self.assertTrue(np.allclose(dac_data, input_data))


class FiberModelTest(unittest.TestCase):
    """Test if the fiber model behaves as expected."""

    def test_fiber(self):
        input_data = read_octave_file("test/data/fiber_input.csv")
        output_data = read_octave_file("test/data/fiber_output.csv")

        fib = models.Fiber()
        fib_data = fib.run(input_data)

        self.assertTrue(np.allclose(fib_data, output_data))


class IQModemModelTest(unittest.TestCase):
    """Test if the iqmodem model behaves as expected."""

    def test_iqmodem_ideal(self):
        input_data = read_octave_file("test/data/iqmodem_input.csv")
        phasor_data = read_octave_file("test/data/iqmodem_phasor.csv")
        output_data = read_octave_file("test/data/iqmodem_output.csv")

        iq = models.IQModem()
        iq.iqi_coef = 1.5
        iq.iqi_filter = 3
        iq.iqi_delay_imbalance = 1.2
        iq.dc_offset = 0.6
        iq_data = iq.run(input_data, phasor_data)

        self.assertTrue(np.allclose(iq_data, output_data))

    def test_iqmodem_filter(self):
        input_data = read_octave_file("test/data/iqmodem_input_filter.csv")
        phasor_data = read_octave_file("test/data/iqmodem_phasor_filter.csv")
        output_data = read_octave_file("test/data/iqmodem_output_filter.csv")

        iq = models.IQModem()
        iq.mode = 'filter'
        iq.iqi_coef = 2.1
        iq.iqi_filter = 1
        iq.iqi_delay_imbalance = 1.8
        iq.dc_offset = 0.2
        iq_data = iq.run(input_data, phasor_data)

        self.assertTrue(np.allclose(iq_data, output_data))

    def test_iqmodem_filter2(self):
        input_data = read_octave_file("test/data/iqmodem_input_filter2.csv")
        phasor_data = read_octave_file("test/data/iqmodem_phasor_filter2.csv")
        output_data = read_octave_file("test/data/iqmodem_output_filter2.csv")

        iq = models.IQModem()
        iq.mode = 'filter'
        iq.iqi_coef = 2.1
        iq.iqi_filter = 7
        iq.iqi_delay_imbalance = 1.8
        iq.dc_offset = 0.2
        iq_data = iq.run(input_data, phasor_data)

        self.assertTrue(np.allclose(iq_data, output_data))

    def test_iqmodem_static(self):
        input_data = read_octave_file("test/data/iqmodem_input_static.csv")
        phasor_data = read_octave_file("test/data/iqmodem_phasor_static.csv")
        output_data = read_octave_file("test/data/iqmodem_output_static.csv")

        iq = models.IQModem()
        iq.mode = 'static'
        iq.iqi_coef = 1.9
        iq.iqi_filter = 3
        iq.iqi_delay_imbalance = 0.6
        iq.dc_offset = 2.4
        iq_data = iq.run(input_data, phasor_data)

        self.assertTrue(np.allclose(iq_data, output_data))


class OscillatorModelTest(unittest.TestCase):
    """Test if the oscillator model behaves as expected."""

    def test_oscillator_ideal(self):
        output_data = [1] * 200

        osc = models.Oscillator()
        osc.mode = 'ideal'
        osc_data = osc.run(200)

        self.assertTrue(np.allclose(osc_data, output_data))

    def test_oscillator_cfo(self):
        numpy.random.seed(45612)
        output_data = read_octave_file("test/data/oscillator_output_cfo.csv").T

        osc = models.Oscillator()
        osc.mode = 'cfo'
        osc.cfo = 200e3
        osc.cfo_std = 320e3
        osc_data = osc.run(200)

        self.assertTrue(np.allclose(osc_data, output_data))

    def test_oscillator_model(self):
        numpy.random.seed(45612)
        output_data = read_octave_file(
            "test/data/oscillator_output_model.csv").T

        osc = models.Oscillator()
        osc.mode = 'model'
        osc.f3db = 250e3
        osc_data = osc.run(200)

        self.assertTrue(np.allclose(osc_data, output_data))

    def test_oscillator_spectrum(self):
        numpy.random.seed(45612)
        output_data = read_octave_file(
            "test/data/oscillator_output_spectrum.csv").T

        osc = models.Oscillator()
        osc.mode = 'spectrum'
        osc.freq = np.linspace(1, 100e6, 10000)
        with open("test/data/spec_response.txt", 'r') as f:
            s = f.readline()
        osc.spec = [float(x) for x in s.split()]
        osc_data = osc.run(200)

        self.assertTrue(np.allclose(osc_data, output_data, rtol=0.05))

    def test_oscillator_variance_spectrum(self):
        numpy.random.seed(45612)

        osc = models.Oscillator()
        osc.mode = 'spectrum'
        osc.freq = np.linspace(1, 100e6, 10000)
        with open("test/data/spec_response.txt", 'r') as f:
            s = f.readline()
        osc.spec = [float(x) for x in s.split()]
        osc_data = osc.variance_phase_spectrum()

        self.assertTrue(np.allclose(osc_data, 2526269123350.2153))


class RadioStripModelTest(unittest.TestCase):
    """Test if the radio strip model behaves as expected."""

    def test_radio_stripe(self):
        numpy.random.seed(45612)

        input_data = read_octave_file("test/data/radiostripe_input.csv")
        output_data = read_octave_file("test/data/radiostripe_output.csv")

        rs = models.RadioStripe()
        rs_data = rs.run(input_data)

        self.assertTrue(np.allclose(rs_data, output_data))

    def test_calibration(self):
        numpy.random.seed(45612)

        input_data = read_octave_file("test/data/radiostripe_input.csv")
        output_data = read_octave_file("test/data/radiostripe_cal_output.csv")

        rs = models.RadioStripe()
        rs.calibrate(input_data, 20)
        rs_data = rs.run(input_data)

        self.assertTrue(np.allclose(rs_data, output_data))


class UtilsTest(unittest.TestCase):
    """Test the utility functions."""

    def test_randconst_default(self):
        """Test randconst for the default QAM16 constellation."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output.csv")

        out, _ = models.utils.randconst(200, 1)

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_qam2(self):
        """Test randconst for QAM2."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_qam2.csv")

        out, _ = models.utils.randconst(200, 1, m=2)

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_qam8(self):
        """Test randconst for QAM8."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_qam8.csv")

        out, _ = models.utils.randconst(200, 1, m=8)

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_qam32(self):
        """Test randconst for QAM32."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_qam32.csv")

        out, _ = models.utils.randconst(200, 1, m=32)

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_qam128(self):
        """Test randconst for QAM128."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_qam128.csv")

        out, _ = models.utils.randconst(200, 1, m=128)

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_qam512(self):
        """Test randconst for QAM512."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_qam512.csv")

        out, _ = models.utils.randconst(200, 1, m=512)

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_psk(self):
        """Test randconst for PSK."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_psk.csv")

        out, _ = models.utils.randconst(200, 1, type='PSK')

        self.assertTrue(np.allclose(out, output_data))

    def test_randconst_spiral(self):
        """Test randconst for SPIRAL."""
        models.utils.global_seed = 45612
        output_data = read_octave_file("test/data/randconst_output_SPIRAL.csv")

        out, _ = models.utils.randconst(200, 1, type='SPIRAL')

        self.assertTrue(np.allclose(out, output_data))

    def test_rrc(self):
        input_data = read_octave_file("test/data/rrc_input.csv")
        output_data = read_octave_file("test/data/rrc_output.csv")

        out = models.utils.rrc(input_data, 0.5)

        self.assertTrue(np.allclose(out, output_data))

    def test_pulseshape_int(self):
        """Test the pulseshape function with an oversampling value which is integer."""
        input_data = read_octave_file("test/data/randconst_output.csv")
        output_data = read_octave_file("test/data/pulseshape_output.csv")

        out = models.utils.pulseshape(input_data, 5, 0.1)

        self.assertTrue(np.allclose(out, output_data))

    def test_pulseshape_float(self):
        """Test the pulseshape function with a floating point oversampling value."""
        input_data = read_octave_file("test/data/randconst_output.csv")
        output_data = read_octave_file("test/data/pulseshape_output2.csv")

        out = models.utils.pulseshape(input_data, 2.5, 0.2)

        self.assertTrue(np.allclose(out, output_data))
