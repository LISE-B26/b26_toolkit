from unittest import TestCase

from b26_toolkit.b26_toolkit.instruments import dcservo_thorlabs_kinesis as DC


class TestTDC001(TestCase):
    def setUp(self):
        self.servo = DC.TDC001()

    def tearDown(self):
        self.servo.__del__()

    def test_connect(self):
        self.assertEqual(self.servo.is_connected, True)

    def test_move(self):
        self.servo.goto_home()
        self.servo.update({'position': 3})
        self.assertAlmostEqual(self.servo.position, 3, places=0)

    def test_vel(self):
        self.servo.update({'velocity': 2.5})
        self.assertEqual(self.servo.velocity, 2.5)