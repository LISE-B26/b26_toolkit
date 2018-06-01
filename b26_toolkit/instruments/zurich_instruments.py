"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""

from pylabcontrol.core import Instrument,Parameter
import numpy as np


# =============== ZURCIH INSTRUMENTS =======================
# ==========================================================
class ZIHF2(Instrument):
    """
instrument class to control the Zurich instrument lockin-amplifier

for instructions about the API to implement more functionality check out Part II Chapter 4 in the manual:
https://www.zhinst.com/sites/default/files/LabOneProgrammingManual_34390_1.pdf
    """

    try:
        import zhinst.utils as utils
        _lib_detected = True
    except ImportError:
        # make a fake ZI instrument
        #import warnings
        #warnings.warn('Unable to import ZI library. Check that it has been installed from the ZI website.')
        _lib_detected = False
    except:
        raise

    _DEFAULT_SETTINGS = Parameter([
        Parameter('freq', 1e6, float, 'frequency of output channel'),
        Parameter('amp', 0.1, float, 'amplitude of output channel (V)'),
        Parameter('sigins',
                  [
                      Parameter('channel', 0, [0,1], 'signal input channel'),
                      Parameter('imp50',False, bool, '50Ohm impedance on (1) or off (0)'),
                      Parameter('ac', False, bool, 'ac coupling on (1) or off (0)'),
                      Parameter('range', 10.0, [0.01, 0.1, 1.0, 10.0], 'range of signal input'),
                      Parameter('diff',  False, bool, 'differential signal on (1) or off (0)')
                   ]
                  ),
        Parameter('sigouts',
                  [
                      Parameter('channel', 0, [0,1], 'signal output channel'),
                      Parameter('on',  False, bool, 'output on (1) or off (0)'),
                      Parameter('add',  False, bool, 'add aux signal on (1) or off (0)'),
                      Parameter('range', 10.0, [0.01, 0.1, 1.0, 10.0], 'range of signal output')
                   ]
                  ),
        Parameter('demods',
                  [
                      Parameter('channel', 0, [0,1], 'demodulation channel'),
                      Parameter('order', 4, int, 'filter order'),
                      Parameter('rate', 10e3, [10e3], 'rate'),
                      Parameter('harmonic', 1, int, 'harmonic at which to demodulate'),
                      Parameter('phaseshift', 0.0, float, 'phaseshift of demodulation'),
                      Parameter('oscselect', 0, [0,1], 'oscillator for demodulation'),
                      Parameter('adcselect', 0, int, 'adcselect not sure what that is!?')
                   ]
                  ),
        Parameter('aux',
                  [
                      Parameter('channel', 0, [0,1], 'auxilary channel'),
                      Parameter('offset', 1.0, float, 'offset in volts')
                   ]
                  )
    ])

    _PROBES = {
            'input1': 'this is the input from channel 1',
            'R': 'the amplitude of the demodulation signal',
            'X': 'the X-quadrature of the demodulation signal',
            'Y': 'the Y-quadrature of the demodulation signal',
            'freq': 'the frequency of the output channel'
        }

    '''
    instrument class to talk to Zurich instrument HF2 lock in ampifier
    '''
    def __init__(self, name = None, settings = None):
        #COMMENT_ME

        self.daq = self.utils.autoConnect(8005,1) # connect to ZI, 8005 is the port number
        self.device = self.utils.autoDetect(self.daq)
        self.options = self.daq.getByte('/%s/features/options' % self.device)

        super(ZIHF2, self).__init__(name, settings)
        # apply all settings to instrument
        self.update(self.settings)

        if self._lib_detected:
            self._is_connected = True


    # ========================================================================================
    # ======= overwrite old_functions from instrument superclass =================================
    # ========================================================================================

    def update(self, settings):
        '''
        updates the internal dictionary and sends changed values to instrument
        Args:
            commands: parameters to be set
        '''
        # call the update_parameter_list to update the parameter list
        super(ZIHF2, self).update(settings)


        def commands_from_settings(settings):
            '''
            converts dictionary to list of  setting, which can then be passed to the zi controler
            :param settings = dictionary that contains the commands
            :return: commands = list of commands, which can then be passed to the zi controler
            '''
            # create list that is passed to the ZI controler


            commands = []



            for key, element in sorted(settings.items()):
                if isinstance(element, dict) and key in ['sigins', 'sigouts', 'demods']:
                    if 'channel' in element:
                        channel = element['channel']
                    else:
                        channel = self.settings[key]['channel']
                    for sub_key, val in sorted(element.items()):
                        if not sub_key == 'channel':
                            commands.append(['/%s/%s/%d/%s'%(self.device, key, channel, sub_key), val])
                        # is the range has been changed we have to reset the amplitude
                        if key == 'sigouts' and sub_key == 'range':
                            channel = self.settings['sigouts']['channel']
                            out_range = self.settings['sigouts']['range']
                            scaled_amplitude = self.settings['amp'] / out_range
                            # I am not sure what the last number means but it is 6 for channel 0 and 7 for channel 1
                            commands.append(['/%s/sigouts/%d/amplitudes/%d' % (self.device, channel, channel + 6),scaled_amplitude])
                elif isinstance(element, dict) and key in ['aux']:
                    if 'channel' in element:
                        channel = element['channel']
                    else:
                        channel = self.settings['aux']['channel']
                    if 'offset' in element:
                        offset = element['offset']
                    else:
                        offset = self.settings['aux']['offset']
                        print(('offset', offset))
                    commands.append(['/%s/AUXOUTS/%d/OFFSET'% (self.device, channel),offset ])
                elif key in ['freq']:
                    channel = self.settings['sigouts']['channel']
                    # commands.append(['/%s/oscs/%d/freq' % (self.device, channel), settings['freq']])
                    commands.append(['/%s/oscs/%d/freq' % (self.device, channel), element])
                elif key in ['amp']:
                    channel = self.settings['sigouts']['channel']
                    out_range = self.settings['sigouts']['range']
                    amplitude = element
                    if amplitude > out_range:
                        # adjust range
                        out_range = [x for x in self._DEFAULT_SETTINGS.valid_values['sigouts']['range'] if amplitude < x][0]
                    # amplitudes are given relative to the range
                    scaled_amplitude = amplitude/out_range
                    commands.append(['/%s/sigouts/%d/range' % (self.device, channel), out_range])
                    # I am not sure what the last number means but it is 6 for channel 0 and 7 for channel 1
                    commands.append(['/%s/sigouts/%d/amplitudes/%d' % (self.device, channel, channel+6), scaled_amplitude])
                elif isinstance(element, dict) == False:
                    commands.append([key, element])


            return commands


        # now we actually apply these newsettings to the hardware
        commands = commands_from_settings(settings)
        if self.is_connected:
            try:
                self.daq.set(commands)
            except RuntimeError:
                print(('runtime error. commands\n{:s}'.format(commands)))
        else:
            print('hardware is not connected, the command to be send is:')
            print(commands)

    def read_probes(self, key):
        '''

        requests value from the instrument and returns it
        Args:
            key: name of requested value

        Returns: reads values from instrument

        '''
        assert key in list(self._PROBES.keys()), "key assertion failed {:s}".format(str(key))

        if key.upper() in ['X', 'Y', 'R']:
            # these values we actually request from the instrument
            data = self.poll(key)
            data = data[key]

        elif key in ['freq']:
            # these values just look up in the parameter settings
            data = self.settings['freq']

        return data

    @property
    def is_connected(self):
        '''
        check if instrument is active and connected and return True in that case
        :return: bool
        '''
        return self._is_connected


    # Poll the value of input 1 for polltime seconds and return the magnitude of the average data. Timeout is in milisecond.
    def poll(self,  variable = 'R', demod_c = 0, pollTime = 0.1, timeout = 500):
        """

        Args:
            variable: string or list of strings, which varibale to poll ('R', 'x', 'y')
            demod_c:
            pollTime: integration time
            timeout: 0.1s could be varibale in the future

        Returns: requested value from instrument as a dictionary {varible: array of values}

        """

        # print('warning! polling from ZI, pollTime and timeout still hardcoded')
        valid_variables = ['R','X','Y']

        if self.is_connected:
            path = '/%s/demods/%d/sample' % (self.device, demod_c)
            self.daq.subscribe(path)
            flat_dictionary_key = True
            data_poll = self.daq.poll(pollTime,timeout,1,flat_dictionary_key)
            data = {}
            for var in ['X','Y']:
                data.update({var: data_poll[path][var.lower()]})
            data.update({'R': np.sqrt(np.square(data['X'])+np.square(data['Y']))})

            if isinstance(variable,str):
                variable = [variable]

            # poll gives an array of values read from the instrument, here we are only instersted in the mean
            return_variable = {k: np.mean(data[k.upper()]) for k in variable}
        else:
            return_variable = None

        return return_variable
