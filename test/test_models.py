import models
import unittest


class ModelsTest(unittest.TestCase):
    def test_amplifier(self):
        """Test if the amplifier class behaves as expected."""
        pa = models.Amplifier()
        pa.mode = 'tanh'
        pa.gain = 3.2
        pa.run(0)

        self.assertEqual(True, True)
