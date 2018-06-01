"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    pylabcontrol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pylabcontrol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.
"""

from pylabcontrol.core import Instrument, Parameter
import visa
import numpy as np
import time

import matplotlib.pyplot as plt


class Oscilloscope(Instrument):
    """
    This class provides a python implementation of the Keysight DSO1024A oscilloscope.
    """

    _INSTRUMENT_IDENTIFIER = 'Agilent Technologies,DSO1024A,CN56301388,00.04.06 SP06'
    # String returned by spectrum analyzer upon querying it with '*IDN?'

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_resource', 'USB0::0x0957::0x0588::CN56301388::0::INSTR', (str),
                      'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
        Parameter('timebase',[
            Parameter('format', 'yt', ['yt', 'xy', 'roll'], 'time base format')
        ]),
        # Parameter('frequency_step', 10e6, float, 'frequency interval of spectrum analyzer frequency range'),
        Parameter('acquisition',[
                      Parameter('type','normal',['normal'], 'acquisition type'),
                      Parameter('count',3, int, 'acquisition count')
                  ]),
        Parameter('waveform', [
            Parameter('mode', 'raw', ['raw', 'normal'], 'waveform mode (not that normal allows only to acquire up to 600 datapoints! raw up tp 10240)'),
            Parameter('points', 10240, int, 'waveform length, max is 10240 (mode raw) or 600 (mode normal)'),
            Parameter('format', 'word', ['word', 'ascii', 'byte'], 'waveform format '),
            Parameter('channel', 1, [1, 2, 3, 4], 'channel from which to read the data'),
            Parameter('timebase', 1, [1,2,5, 10, 20, 50, 100, 200, 500], 'timebase: units per devision'),
            Parameter('timebase_unit', 'ms', ['ns', 'us','ms', 's'], 'timebase: units per devision'),
            Parameter('vert_scale', '5E-2', ['2E-3', '5E-3', '1E-2', '2E-2', '5E-2', '1E-1', '2E-1', '5E-1'], 'vert_scale volts per division (2 mV to 10 V)'),
            #Parameter('offset', 0.0, float, 'vertical scale off set in V'),
            Parameter('offset', '-1E-2', ['-5E-1', '-2E-1', '-1E-1', '-5E-2', '-2E-2', '-1E-2', '-5E-3', '-2E-3', '2E-3', '5E-3', '1E-2', '2E-2', '5E-2', '1E-1', '2E-1', '5E-1'],
                      'offset on channel in volts'),
            Parameter('probe', '1X', ['1X', '10X'], 'multipliation factor for y scale')
        ]),
        Parameter('trigger', [
            Parameter('on', True, bool, 'trigger on or off'),
            Parameter('channel', 2, [1, 2, 3, 4], 'channel from which to trigger'),
            Parameter('mode', 'edge_pos', ['edge_pos', 'edge_neg'], 'trigger mode (psotive edge, negative edge)'),
            Parameter('level', 0.0, [0.0, 0.1, 1], 'trigger level (I think in fractions of a devision between -6 and 6 could also be in volt)'),

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

        # keep track of when the instrument was updated last to prevent sending requests to frequently
        self._last_update_time = time.time()

        rm = visa.ResourceManager()

        # todo: JG 20170623 implement proper error handling when insturment is not connected.
        self.osci = rm.open_resource(self.settings['visa_resource'])
        self.osci.read_termination = '\n'
        self.osci.timeout = self.settings['connection_timeout']
        self.osci.inputbuffersize = 1000
        self.osci.ByteOrder = 'littleEndian'

        self.osci.write('*RST') #Places the oscilloscope in the factory default setup state.
        self._wait_for_osci()
        self.update(self.settings)
        # except:
        #     raise

    def reset(self):
        """
        resets the oscilloscope to the default values
        Returns:

        """
        self.osci.write('*RST') #Places the oscilloscope in the factory default setup state.

    def update(self, settings):
        """
        updates the instrument parameters

        Args:
            settings: dictionary that contains the parameter indentifiers (keys) and the new parameters values (value)

        """
        super(Oscilloscope, self).update(settings)

        if 'timebase' in settings:
            # self._wait_for_osci()
            # self._set_timebase(settings['timebase'])
            self.osci.write(':TIMEBASE:MODE MAIN')
            self.osci.write(':TIMEBASE:FORM ' + settings['timebase']['format'].upper())

        if 'acquisition' in settings:
            # self._wait_for_osci()
            # self._set_acquisition(settings['acquisition'])
            self.osci.write(':ACQUIRE:TYPE ' + settings['acquisition']['type'].upper())
            self.osci.write(':ACQUIRE:COUNT ' + str(settings['acquisition']['count']))

        if 'waveform' in settings:
            channel = str(settings['waveform']['channel'])
            waveform = settings['waveform']
            # self._wait_for_osci()
            # self._set_waveform(settings['waveform'])
            self.osci.write(':WAV:POINTS:MODE ' + waveform['mode'].upper())
            self.osci.write(':WAV:POINTS ' + str(waveform['points']))
            self.osci.write(':WAV:FORMAT ' + waveform['format'].upper())
            self.osci.write(':WAV:SOURCE CHAN' + channel)
            self.osci.write(':TIM:SCAL ' + self.time_base_to_nr3(waveform['timebase'], waveform['timebase_unit']))
            self.osci.write(':CHAN' + channel + ':PROB ' + waveform['probe']) #PROBE MUST BE SET BEFORE SCALE
            self.osci.write(':CHAN' + channel + ':SCAL ' + waveform['vert_scale'])
            self.osci.write(':CHAN' + channel + ':OFFS ' + waveform['offset'])#'{:0.02e}'.format(waveform['offset']))

        if 'trigger' in settings:
            channel = str(settings['trigger']['channel'])
            self.osci.write(':TRIG:PULS:SOUR ' + channel)
            if 'mode' in settings['trigger']:
                if settings['trigger'] is 'edge_pos':
                    self.osci.write(':TRIG:EDGE:SLOP POS')
                if settings['trigger'] is 'edge_neg':
                    self.osci.write(':TRIG:EDGE:SLOP NEG')

            self.osci.write(':TRIG:EDGE:LEV ' + str(settings['trigger']['level']))


    def acq_time(self):
        """
        estimates the acquisition time
        Returns:

        """
        waveform = self.settings['waveform']
        total_time = float(self.time_base_to_nr3(waveform['timebase'], waveform['timebase_unit']))*10
        return total_time


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

    def get_timetrace(self):
        """
        reads a time trace from the intruments
        Returns:

        """
        self.osci.write(':SINGLE') # start a single acquisition
        self.osci.write(':FORCetrig') # force trigger to make sure that the acquisition starts


        # estimate acquisition time and wait until acquisition is finished
        # JG: This seems to work, just wait enough time so that the oscilloscope can acquire the data
        # output counter to terminal to show that the program didn't freeze
        total_time = self.acq_time()
        total_time = int(total_time )*10 # multiply by 10 because the wait loop time is 100ms
        total_time = int(total_time * 1.2+1) # give another percent margin, this is empirical

        for i in range(total_time):
            print(('waiting {:d}/{:d}'.format(i, total_time)))
            time.sleep(0.1)

        operationComplete = bool(self.osci.query('*OPC?'))
        print(('operationComplete', operationComplete))


        # Get the preamble block
        preambleBlock = self.osci.query(':WAV:PREAMBLE?')
        # preable contains the curren settings
        #   FORMAT        : int16 - 0 = WORD, 1 = BYTE, 2 = ASCII.
        #   TYPE          : int16 - 0 = NORMAL, 1 = PEAK DETECT, 2 = AVERAGE
        #   POINTS        : int32 - number of data points transferred.
        #   COUNT         : int32 - 1 and is always 1.
        #   XINCREMENT    : float64 - time difference between data points.
        #   XORIGIN       : float64 - always the first data point in memory.
        #   XREFERENCE    : int32 - specifies the data point associated with x-origin.
        #   YINCREMENT    : float32 - voltage diff between data points.
        #   YORIGIN       : float32 - value is the voltage at center screen.
        #   YREFERENCE    : int32 - specifies the data point where y-origin occurs

        # convert into dictionary
        preamble = {k:float(v) for k, v in zip(['format', 'type', 'points', 'count', 'xincrement', 'xorigin', 'xreference', 'yincrement', 'yorigin', 'yreference'], preambleBlock.split(','))}



        format = str(self.settings['waveform']['format']).lower()
        # depending on the setting we get differnet data back
        if format == 'ascii':
            print('datatype ascii')
            print('JG: WARNING ascii NOT TESTED!')

            # send command to read data
            # raw_data = self.osci.query(':WAV:DATA?')
            raw_data = self.osci.query_ascii_values(':WAV:DATA?')
            # the first charactars are some metadata
            prefix = raw_data[:11]
            data = raw_data[11:].split(',')
            # data = [float(d) for d in data]
            data = np.array(data)

        elif format == 'word':
            print('datatype word')
            # try until we get the data
            max_attempts = 100
            count = 0
            raw_data = []
            while len(raw_data) == 0 and count<max_attempts:
                raw_data = self.osci.query_binary_values(':WAV:DATA?', datatype='H')
                print(('data empty retry, attempt {:02d}/{:02d}'.format(count, max_attempts)))
                count += 1
                time.sleep(0.1)
            data = raw_data

        elif format == 'byte':
            print('datatype byte')
            print('JG: WARNING BYTE NOT TESTED!')
            raw_data = self.osci.query_binary_values(':WAV:DATA?')
            data = raw_data
        else:
            print('WARNING UNKNOWN DATA FORMAT')

        if len(data)>0:
            dt = float(self.time_base_to_nr3(self.settings['waveform']['timebase'], self.settings['waveform']['timebase_unit']))*10/len(data)
        else:
            dt = 1
        # add more meta data to the preamble
        preamble['dt']= dt
        preamble['vert_scale'] =float(self.settings['waveform']['vert_scale'])*float(self.settings['waveform']['probe'].split('X')[0]) #  vertical scale of osci
        return data, preamble


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


    @staticmethod
    def time_base_to_nr3(n, s):
        """
        converts the time base into the nr3 format
        Args:
            n: time base value (1, 2, 5)
            s: time base (ns, us, ms, s)

        Returns:
            a string in nr3 format that represents the time for example, 1.0E-9, 2.0E-9, 5.0E-9, ... 1.0E+00, 2.0E+00, 5.0E+00)

        """

        if s == 'ns':
            nr3 = 'E-9'
        if s == 'us':
            nr3 = 'E-6'
        if s == 'ms':
            nr3 = 'E-3'
        if s == 's':
            nr3 = 'E+00'

        return str(n) + '.0' + nr3

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
        print(oscil.is_connected())
        # print oscil.mode
        # oscil.mode = 'TrackingGenerator'
        # print oscil.mode

        print('=============')

        print((oscil.settings))

        data, preambleBlock = oscil.get_timetrace()

        dt = preambleBlock['dt']
        time = dt*np.arange(len(data))
        print(('data', data))
        plt.plot(time, data, '-x')

        plt.show()
