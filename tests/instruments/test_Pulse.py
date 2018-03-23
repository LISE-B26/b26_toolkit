from unittest import TestCase

import numpy as np

from b26_toolkit.src.instruments import  Pulse
from collections import namedtuple

PulseTuple = namedtuple('Pulse', ('channel_id', 'start_time', 'duration'))

class TestPulse(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_compare_to_tuple(self):

        start = 0
        duration = 10
        pt = PulseTuple('p1', start, duration)

        p = Pulse('p1', start, duration)

        assert p.channel_id == pt.channel_id
        assert p.start_time == pt.start_time
        assert p.duration == pt.duration


    def test_overlapping(self):


        p1 = Pulse('p1', 1, end_time=3)
        p2 = Pulse('p2', 4, duration=2)

        assert Pulse.is_overlapping(p1, p2, dead_time=0) is False
        assert Pulse.is_overlapping(p1, p2, dead_time=1) is False
        assert Pulse.is_overlapping(p1, p2, dead_time=2) is True
