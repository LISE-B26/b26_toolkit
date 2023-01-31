
from pylabcontrol.core import Parameter, Script
from b26_toolkit.instruments import AFG3022C_02, ANC300


class LoadInstrument(Script):
    """
    Import this script to load some instrument.
    """

    _DEFAULT_SETTINGS = [
        Parameter('this', 'script', ['script', 'was simply imported', 'to load some instrument'], 'change the instrument to whatever you want in the code')
    ]
    _SCRIPTS = {}
    _INSTRUMENTS = {'ANC300': ANC300}
    #_INSTRUMENTS = {'afg_2': AFG3022C_02}
