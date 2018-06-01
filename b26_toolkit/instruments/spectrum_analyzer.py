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

from pylabcontrol.core import Instrument, Parameter
import visa
import numpy as np
import time


class SpectrumAnalyzer(Instrument):
    """
    This class provides a python implementation of the Keysight N9320B 9kHz-3.0GHz spectrum analyzer
    with trigger generator.
    """

    _INSTRUMENT_IDENTIFIER = 'Keysight Technologies,N9320B,CN0323B356,0B.03.58'
    # String returned by spectrum analyzer upon querying it with '*IDN?'

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_resource', 'USB0::0x0957::0xFFEF::CN0323B356::INSTR', (str),
                      'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
        Parameter('start_frequency', 9e3, float, 'start frequency of spectrum analyzer frequency range'),
        Parameter('mode', 'SpectrumAnalyzer', ['SpectrumAnalyzer', 'TrackingGenerator'],
                  'switches between normal spectrum analyzer mode or spectrum analyzer PLUS output, '
                  'i.e., tracking generator'),
        Parameter('stop_frequency', 3e9, float, 'stop frequency of spectrum analyzer frequency range'),
        # Parameter('frequency_step', 10e6, float, 'frequency interval of spectrum analyzer frequency range'),
        Parameter('output_on', False, bool, 'toggles the tracking generator'),
        Parameter('connection_timeout', 1000, int, 'the time to wait for a response '
                                                   'from the spectrum analyzer with each query (units??)'),
        Parameter('output_power', -20.0, float, 'the output power (in dBm) of the tracking generator')
        ])

    _PROBES = {'start_frequency': 'the lower bound of the frequency sweep',
               'stop_frequency': 'the upper bound of the frequency sweep',
               'trace': 'the frequency sweep of the inputted signal',
               'tracking_generator': 'checks if the tracking generator is on',
               'bandwidth': 'the curent bandwidth of the spectrum analyzer',
               'output_power': 'the power of the tracking generator',
               'mode': 'Spectrum Analyzer Mode or Tracking Generator Mode'}

    def __init__(self, name='SpectrumAnalyzer', settings={}):
        """

        Args:
            name (str): optional name of instance of class
            settings (list): list of other values to initialize class with

        """

        super(SpectrumAnalyzer, self).__init__(name, settings)

        self._last_update_time = time.time()

        rm = visa.ResourceManager()

        # todo: JG 20170623 implement proper error handling when insturment is not connected.
        # try:
        self.spec_anal = rm.open_resource(self.settings['visa_resource'])
        self.spec_anal.read_termination = '\n'
        self.spec_anal.timeout = self.settings['connection_timeout']
        self.spec_anal.write('*RST') #Places the oscilloscope in the factory default setup state.
        self._wait_for_spec_anal()
        self.update({'mode':'SpectrumAnalyzer'})
        # except:
        #     raise

    def update(self, settings):
        """
        updates the instrument parameters

        Args:
            settings: dictionary that contains the parameter indentifiers (keys) and the new parameters values (value)

        """
        super(SpectrumAnalyzer, self).update(settings)

        # set mode first
        if 'mode' in settings:
            self._wait_for_spec_anal()
            self._set_mode(settings['mode'])
            print(('mode',settings['mode']))
            # since changes in the output_power are not applied to the instruments when in SpectrumAnalyzer mode, we make sure that is is updated once switched back to TrackingAnalyzer
            if settings['mode'] == 'TrackingGenerator':
                self.update({'output_power':self.settings['output_power']})

        if 'start_frequency' in settings:
            assert 9e3 <= settings[
                'start_frequency'] <= 3e9, "start frequency must be between 0 and 3e9, you tried to set it to {0}!".format(
                settings['start_frequency'])
            self._wait_for_spec_anal()
            self._set_start_frequency(settings['start_frequency'])

        # if 'frequency_step' in settings:
        #     self._wait_for_osci()
        #     self._set_frequency_step(settings['frequency_step'])

        if 'stop_frequency' in settings:
            assert 9e3 <= settings[
                'stop_frequency'] <= 3e9, "start frequency must be between 0 and 3e9, you tried to set it to {0}!".format(
                settings['stop_frequency'])
            self._wait_for_spec_anal()
            self._set_stop_frequency(settings['stop_frequency'])

        if 'output_on' in settings:
            self._wait_for_spec_anal()
            self._toggle_output(settings['output_on'])

            # todo: JG 201700623 check if the power is set properly (see comment in "if 'mode' in settings:" (~ line 87)
            print('warning! output turned on.\ncheck if the output power is corresponds to the value in the settings!\nCheck code for more details' )


        if 'output_power' in settings and self.settings['output_on']:
            self._wait_for_spec_anal()
            self._set_output_power(settings['output_power'])

        # for key, value in settings.iteritems():
            # if key == 'start_frequency':
            #     assert 0.0 < value < 3e9, \
            #         "start frequency must be between 0 and 3e9, you tried to set it to {0}!".format(value)
            #     self._set_start_frequency(value)
            # elif key == 'stop_frequency':
            #     assert 0.0 < value < 3e9, \
            #         "stop frequency must be between 0 and 3e9, you tried to set it to {0}!".format(value)
            #     self._set_stop_frequency(value)
            # elif key == 'output_on':
            #     self._toggle_output(value)
            # elif key == 'output_power':
            #     self._set_output_power(value)
            # elif key == 'mode':
            #     self._set_mode(value)
            # else:
            #     message = '{0} is not a parameter of {1}'.format(key, self.name)

    def read_probes(self, probe_name):
        self._wait_for_spec_anal()

        if probe_name == 'start_frequency':
            return self._get_start_frequency()
        elif probe_name == 'stop_frequency':
            return self._get_stop_frequency()
        elif probe_name == 'trace':
            return self._get_trace()
        elif probe_name == 'output_on':
            return self._is_output_on()
        elif probe_name == 'bandwidth':
            return self._get_bandwidth()
        elif probe_name == 'output_power':
            return self._get_output_power()
        elif probe_name == 'mode':
            return self._get_mode()
        else:
            message = 'no probe with that name exists!'
            raise AttributeError(message)

    def is_connected(self):
        """
        Checks if the instrument is connected.
        Returns: True if connected, False otherwise.

        """
        identification = self.spec_anal.query('*IDN?')
        return identification == self._INSTRUMENT_IDENTIFIER

    def _set_start_frequency(self, start_freq):
        #COMMENT_ME
        self.spec_anal.write('SENS:FREQ:START ' + str(start_freq))

    def _set_frequency_span(self, frequency_step):
        #COMMENT_ME
        self.spec_anal.write('SENS:FREQ:SPAN ' + str(frequency_step))

    def _get_start_frequency(self):
        #COMMENT_ME
        return float(self.spec_anal.query('SENS:FREQ:START?\n'))

    def _set_stop_frequency(self, stop_freq):
        #COMMENT_ME
        self.spec_anal.write('SENS:FREQ:STOP ' + str(stop_freq))

    def _get_stop_frequency(self):
        #COMMENT_ME
        return float(self.spec_anal.query('SENS:FREQ:STOP?\n'))

    def _toggle_output(self, state):
        #COMMENT_ME

        if state:
            assert self._get_mode() == 'TrackingGenerator', "output can't be on while in SpectrumAnalyzer mode"
            self.spec_anal.write('OUTPUT 1')
        elif not state:
            self.spec_anal.write('OUTPUT 0')

    def _is_output_on(self):
        #COMMENT_ME
        if self.mode == 'SpectrumAnalyzer':
            return False
        elif self.mode == 'TrackingGenerator':
            return bool(int(self.spec_anal.query('OUTPUT:STATE?')))

    def _get_mode(self):
        #COMMENT_ME
        mode_response = str(self.spec_anal.query('CONFIGURE?')).strip()
        if mode_response == 'SAN':
            return 'SpectrumAnalyzer'
        elif mode_response == 'TGEN':
            return 'TrackingGenerator'

    def _set_mode(self, mode):
        #COMMENT_ME
        if mode == 'TrackingGenerator':
            self.spec_anal.write('CONFIGURE:TGENERATOR')
        elif mode == 'SpectrumAnalyzer':
            self.output_on = False
            self.spec_anal.write('CONFIGURE:SANALYZER')

    def _get_trace(self):
        #COMMENT_ME
        amplitudes = [float(i) for i in str(self.spec_anal.query('TRACE:DATA? TRACE1'+ ';*OPC?')).rstrip(';1').split(',')]
        num_points = len(amplitudes)
        frequencies = np.linspace(start=self.start_frequency, stop=self.stop_frequency,
                                  num=num_points).tolist()
        # return [(frequencies[i], amplitudes[i])for i in range(num_points)]

        return {'frequencies':frequencies, 'amplitudes':amplitudes}

    def _get_bandwidth(self):
        #COMMENT_ME
        return float(self.spec_anal.query('BANDWIDTH?'))

    def _get_output_power(self):
        #COMMENT_ME
        return float(self.spec_anal.query('SOURCE:POWER?'))

    def _set_output_power(self, power):
        #COMMENT_ME
        assert self.mode == 'TrackingGenerator', "mode need to be 'TrackingGenerator' to change power"

        return self.spec_anal.write('SOURCE:POWER ' + str(power))

    def __del__(self):
        #COMMENT_ME
        self._wait_for_spec_anal()
        self._set_mode('SpectrumAnalyzer')
        self.spec_anal.close()

    def _wait_for_spec_anal(self):
        #COMMENT_ME

        if self._last_update_time - time.time() < 1.0:
            time.sleep(1)

        self._last_update_time = time.time()


if __name__ == '__main__':

        # sett = {
        #     "settings": {
        #         "output_power": -20.0,
        #         "stop_frequency": 3000000000.0,
        #         "start_frequency": 0.0,
        #         "connection_timeout": 1000,
        #         "mode": "SpectrumAnalyzer",
        #         "output_on": False,
        #         "visa_resource": "USB0::0x0957::0xFFEF::CN0323B356::INSTR"
        #     }
        # }
        #
        print('creat spectrum analyzer instance:')
        spec_anal = SpectrumAnalyzer()
        print(spec_anal.is_connected())
        print(spec_anal.mode)
        spec_anal.mode = 'TrackingGenerator'
        print(spec_anal.mode)

        print('=============')

        print((spec_anal.settings, type(spec_anal.settings)))
