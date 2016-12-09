# clr is python for .net
import clr # run pip install pythonnet
import sys, os
import time
from PyLabControl.src.core.read_write_functions import get_config_value

# dll_path = get_config_value('KINESIS_DLL_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
dll_path = 'C:\Program Files (x86)\Newport\MotionControl\SMC100\Bin'
sys.path.insert(0,dll_path)
from PyLabControl.src.core import Parameter, Instrument
import System

clr.AddReference('Newport.SMC100.CommandInterface')
import CommandInterfaceSMC100

#status codes
MOVING = '28'
DONE_MOVING = '33'

class SMC100(Instrument):
    '''

    '''

    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM9', str, 'serial number written on device'),
        Parameter('position', 25, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('height_lower_limit', 0, float, 'lowest position servo can move to (in mm)')
    ])

    def __init__(self, name = None, settings = None):
        super(SMC100, self).__init__(name, settings)

        self.SMC = CommandInterfaceSMC100.SMC100()
        result = self.SMC.OpenInstrument(self.settings['port'])
        if result == -1:
            print('Failed to open device on port' + str(self.settings['port']))
            raise

    def update(self, settings):
        '''
        Updates internal settings, as well as the position and velocity set on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        '''
        super(SMC100, self).update(settings)
        for key, value in settings.iteritems():
            if key == 'position':
                self._set_position(value)

    @property
    def _PROBES(self):
        return{
            'position': 'servo position in mm'
        }

    def read_probes(self, key):
        '''
        requestes value from the instrument and returns it
        Args:
            key: name of requested value

        Returns: reads values from instrument

        '''
        assert key in self._PROBES.keys()
        assert isinstance(key, str)

        if key == 'position':
            return self._get_position()

    def __del__(self):
        '''
        Cleans up TDC001 connection
        :PostState: TDC001 is disconnected
        '''
        self.SMC.CloseInstrument()

    def _get_position(self):
        i_ref = -1
        s_ref = ''
        result, response, errString = self.SMC.TP(1, i_ref, s_ref)
        if result == -1:
            print('ERROR: ' + errString)
            raise
        return response

    def _set_position(self, pos):
        if pos < self.settings['height_lower_limit']:
            raise ValueError('cannot set position below height_lower_limit')
        s_ref = ''
        result, errString = self.SMC.PA_Set(1, pos, s_ref)
        if result == -1:
            print('ERROR: ' + errString)
            raise
        #block until movement is done
        while True:
            result, ErrorCode, StatusCode, errString = self.SMC.TS(1, s_ref, s_ref,s_ref)
            if result == -1:
                print('ERROR: ' + errString)
                raise
            if StatusCode == DONE_MOVING:
                break
            time.sleep(.1)

# a = SMC100()
# a._set_position(25)
# a._set_position(10)
# print(a._get_position())
# inst, failed = Instrument.load_and_append({'SMC100': SMC100})
# print(inst)
# print(failed)
