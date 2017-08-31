"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    PyLabControl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PyLabControl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PyLabControl.  If not, see <http://www.gnu.org/licenses/>.
"""

from PyLabControl.src.core import Instrument, Parameter
import visa
import numpy as np
import time


class Oscilloscope(Instrument):
    """
    This class provides a python implementation of the Keysight DSO1024A oscilloscope.
    """

    _INSTRUMENT_IDENTIFIER = u'Agilent Technologies,DSO1024A,CN56301388,00.04.06 SP06'
    # String returned by spectrum analyzer upon querying it with '*IDN?'

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_resource', 'USB0::0x0957::0x0588::CN56301388::0::INSTR', (str),
                      'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
        Parameter('channel', 1, [1, 2, 3, 4], 'channel from which to read the data'),
        Parameter('timebase', 'main', ['main'], 'time base mode'),
        # Parameter('frequency_step', 10e6, float, 'frequency interval of spectrum analyzer frequency range'),
        Parameter('acquisition',[
                      Parameter('type','normal',['normal'], 'acquisition type'),
                      Parameter('count', 1, int, 'acquisition count')
                  ]),
        Parameter('waveform', [
            Parameter('mode', 'raw', ['raw'], 'waveform mode'),
            Parameter('points', 5000, int, 'waveform length')
        ]),
        Parameter('connection_timeout', 1000, int, 'the time to wait for a response from the oscilloscope with each query (units??)'),
        ])

    _PROBES = {}
    # _PROBES = {'start_frequency': 'the lower bound of the frequency sweep',
    #            'stop_frequency': 'the upper bound of the frequency sweep',
    #            'trace': 'the frequency sweep of the inputted signal',
    #            'tracking_generator': 'checks if the tracking generator is on',
    #            'bandwidth': 'the curent bandwidth of the spectrum analyzer',
    #            'output_power': 'the power of the tracking generator',
    #            'mode': 'Spectrum Analyzer Mode or Tracking Generator Mode'}

    def __init__(self, name='Oscilloscope', settings={}):
        """

        Args:
            name (str): optional name of instance of class
            settings (list): list of other values to initialize class with

        """

        super(Oscilloscope, self).__init__(name, settings)

        self._last_update_time = time.time()

        rm = visa.ResourceManager()

        # todo: JG 20170623 implement proper error handling when insturment is not connected.
        # try:
        self.osci = rm.open_resource(self.settings['visa_resource'])
        self.osci.read_termination = '\n'
        self.osci.timeout = self.settings['connection_timeout']
        self.osci.write('*RST')
        self._wait_for_osci()
        # self.update({'mode':'SpectrumAnalyzer'})
        # except:
        #     raise

    def update(self, settings):
        """
        updates the instrument parameters

        Args:
            settings: dictionary that contains the parameter indentifiers (keys) and the new parameters values (value)

        """
        super(Oscilloscope, self).update(settings)

        # set mode first
        if 'channel' in settings:
            self._wait_for_osci()
            # self._set_channel(settings['channel'])
            self.osci.write(':WAVEFORM:SOURCE CHAN' + str(settings['channel']))

        if 'timebase' in settings:
            # self._wait_for_osci()
            # self._set_timebase(settings['timebase'])
            self.osci.write(':TIMEBASE:MODE ' + settings['timebase'].upper())

        if 'acquisition' in settings:
            # self._wait_for_osci()
            # self._set_acquisition(settings['acquisition'])
            self.osci.write(':ACQUIRE:TYPE ' + settings['acquisition']['type'].upper())
            self.osci.write(':ACQUIRE:COUNT ' + settings['acquisition']['count'].upper())

        if 'waveform' in settings:
            # self._wait_for_osci()
            # self._set_waveform(settings['waveform'])
            self.osci.write(':WAV:POINTS:MODE ' + settings['waveform']['mode'].upper())
            self.osci.write(':WAV:POINTS: ' + str(settings['waveform']['points']))

    def read_probes(self, probe_name):
        print('NOT IMPLEMENTED')
        pass
        # self._wait_for_osci()
        #
        # if probe_name == 'start_frequency':
        #     return self._get_start_frequency()
        # elif probe_name == 'stop_frequency':
        #     return self._get_stop_frequency()
        # elif probe_name == 'trace':
        #     return self._get_trace()
        # elif probe_name == 'output_on':
        #     return self._is_output_on()
        # elif probe_name == 'bandwidth':
        #     return self._get_bandwidth()
        # elif probe_name == 'output_power':
        #     return self._get_output_power()
        # elif probe_name == 'mode':
        #     return self._get_mode()
        # else:
        #     message = 'no probe with that name exists!'
        #     raise AttributeError(message)

    def is_connected(self):
        """
        Checks if the instrument is connected.
        Returns: True if connected, False otherwise.

        """
        identification = self.osci.query('*IDN?')
        return identification == self._INSTRUMENT_IDENTIFIER


    # def _set_source_channel(self, chan):
    #     """
    #     sets the source channel for acquisition
    #
    #     Returns:
    #
    #     """
    #     self.osci.write(':WAVEFORM:SOURCE CHAN' + str(chan))
    #
    # def _set_timebase(self, timebase):
    #     """
    #     sets the source channel for acquisition
    #
    #     Returns:
    #
    #     """
    #     self.osci.write(':TIMEBASE:MODE ' + timebase.upper())

    def _get_timetrace(self):
        """
        reads a time trace from the intruments
        Returns:

        """

        self.osci.write(':STOP')

        self.osci.write(':DIGITIZE CHAN1') # JG: tmp

        operationComplete = self.osci.query('*OPC?')
        # Get the data back as a WORD (i.e., INT16), other options are ASCII and BYTE
        self.osci.write(':WAVEFORM:FORMAT WORD')
        print('adsasda', operationComplete)
        # Get the preamble block
        preambleBlock = self.osci.query(':WAVEFORM:PREAMBLE?')

        # send command to read data
        raw_data = self.osci.write(':WAV:DATA?')
        print(raw_data)
        # amplitudes = [float(i) for i in str(self.osci.query('TRACE:DATA? TRACE1' + ';*OPC?')).rstrip(';1').split(',')]
        # num_points = len(amplitudes)
        # frequencies = np.linspace(start=self.start_frequency, stop=self.stop_frequency,
        #                           num=num_points).tolist()
        # # return [(frequencies[i], amplitudes[i])for i in range(num_points)]

        # return {'frequencies':frequencies, 'amplitudes':amplitudes}


    def __del__(self):
        #COMMENT_ME
        # self._wait_for_osci()
        # self._set_mode('SpectrumAnalyzer')
        self.osci.close()

    def _wait_for_osci(self):
        #COMMENT_ME
        print('waiting')
        if self._last_update_time - time.time() < 1.0:
            print('sleeping')
            time.sleep(1)

        self._last_update_time = time.time()
        print('-- waited')


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
        print('create oscilloscope instance:')
        oscil = Oscilloscope()
        print oscil.is_connected()
        # print oscil.mode
        # oscil.mode = 'TrackingGenerator'
        # print oscil.mode

        print('=============')

        print(oscil.settings)
