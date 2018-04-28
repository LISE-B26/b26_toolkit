from unittest import TestCase

from pylabcontrol.src.core import Script


class TestScriptDummy(TestCase):


    def setUp(self):
        pass

    def test_xy8_double_init(self):


        updated_scripts, load_failed, updated_instruments = Script.load_and_append({'XY8': 'XY8_double_init'},
                                                                                   package='b26_toolkit')
        xy8 = updated_scripts['XY8']

        xy8.update({'Tracking': {
            'on/off': False}})  # turn off tracking because this will cause an error if we don't run findnv
        print(xy8)

        xy8.is_valid()









