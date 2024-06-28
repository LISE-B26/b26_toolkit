from b26_toolkit.instruments import ArduinoZero
from pylabcontrol.core import Parameter, Script


class LoadInstrument(Script):
    """
    Import this script to load some instrument.
    Actually you don't even need this. You can choose to export an instrument and then load it just like any script.
    The instrument just needs to be in the __init__ file
    """

    _DEFAULT_SETTINGS = [
        Parameter('this', 'script', ['script', 'was simply imported', 'to load some instrument'],
                  'change the instrument to whatever you want in the code')
    ]
    _SCRIPTS = {}
    _INSTRUMENTS = {'arduino': ArduinoZero}

