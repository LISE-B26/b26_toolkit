from unittest import TestCase

from b26_toolkit.b26_toolkit.instruments import CryoStation


class TestCryoStation(TestCase):
    def test_init(self):

        inst = CryoStation()

        print((inst._PROBES))
        print((inst.platform_temp))
        print((inst.stage_1_temp))


