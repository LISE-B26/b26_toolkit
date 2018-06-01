from unittest import TestCase

from b26_toolkit.b26_toolkit.instruments import MaestroBeamBlock, MaestroController


class TestMaestroBeamBlock(TestCase):

    def test_init(self):
        maestro = MaestroController()
        print(maestro)
        test = MaestroBeamBlock(maestro)

        test.is_connected