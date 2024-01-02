
from pylabcontrol.core import Parameter, Script
from b26_toolkit.instruments import AFG3022C_02, ANC300, ArduinoZero, AFG3022C
from pylabcontrol.core import Parameter, Script
from b26_toolkit.instruments import RFGenerator, B26PulseBlaster


class LoadInstrument(Script):
    """
    Import this script to load some instrument.
    """

    _DEFAULT_SETTINGS = [
        Parameter('this', 'script', ['script', 'was simply imported', 'to load some instrument'], 'change the instrument to whatever you want in the code')
    ]
    _SCRIPTS = {}
    #_INSTRUMENTS = {'arduino': ArduinoZero}
    _INSTRUMENTS = {'arduino': ArduinoZero}

