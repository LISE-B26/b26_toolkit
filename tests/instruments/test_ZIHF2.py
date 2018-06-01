from unittest import TestCase

from b26_toolkit.b26_toolkit.instruments import ZIHF2


class TestZIHF2(TestCase):
    def test_init(self):
        test = ZIHF2('my zi')
        test = ZIHF2('my ZI', {'freq':1.})
    def test_update_request(self):

        test = ZIHF2()
        test.update({'sigins':
                                    {'diff': False,
                                     'ac': True,
                                     'imp50': False,
                                     'range': 10,
                                     'channel': 0}
                                })

        self.assertEqual(test.settings['sigins'], {'diff': False,
                                     'ac': True,
                                     'imp50': False,
                                     'range': 10,
                                     'channel': 0}
                         )


        test.update({'freq':  100.})
        self.assertEqual(test.freq, 100.)


        test.settings['freq'] = 101.

        self.assertEqual(test.freq, 101.)
