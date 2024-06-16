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
from struct import unpack
import pyvisa
import numpy as np
import time

import matplotlib.pyplot as plt

# class Oscilloscope(Instrument):

class RigolOscilloscope(Instrument):

    """
    This class provides a python implementation of the Keysight DSO1024A oscilloscope.
    """

    # String returned by spectrum analyzer upon querying it with '*IDN?'
    _INSTRUMENT_IDENTIFIER = 'RIGOL TECHNOLOGIES,MSO5104,MS5A241503457,00.01.03.02.02'

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_resource', 'USB0::0x1AB1::0x0515::MS5A241503457::INSTR', str,
                      'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
        Parameter('timebase_format', 'MAIN', ['MAIN', 'XY', 'ROLL'], 'time base format'),
        Parameter('acq_type', 'NORM', ['NORM', 'AVER', 'PEAK', 'HRES'], 'acquisition type'),
        Parameter('acq_count', 1, int, 'acquisition count'),
        Parameter('acq_memory_depth', '10k', ['AUTO', '1k', '10k','100k','1M','5M','10M','25M','50M','100M','200M'],'number of waveform points per trigger'),
        Parameter('channel', 'CHAN1', ['CHAN1', 'CHAN2', 'CHAN3', 'CHAN4'], 'channel from which to read the data'),
        Parameter('mode', 'RAW', ['NORM', 'MAX', 'RAW'], 'waveform mode (not that normal allows only to acquire up to 1000 datapoints!)'),
        Parameter('format', 'BYTE', ['WORD', 'ASC', 'BYTE'], 'waveform format '),
        Parameter('timebase', 2.0e1, float, 'timebase: units per division'),
        Parameter('offset', 0.0, float, 'offset on channel in volts. range is +/- 1V, +/1 30V, +/- 100V depending on vertical scale'),
        Parameter('probe', 1, int,
                      'multiplication factor for y scale - can only be 1s, 2s, 5s from 1e-4 to 1e4'),
        Parameter('vert_scale', 500.0e-3, float, 'volts per division (500 uV to 10 V)'),
        Parameter('connection_timeout', 9., float, 'the time to wait for a response from the oscilloscope with each query in seconds (16 ns to 10s'),
        ])

    _COMMANDS = {
        'timebase_format': ':TIM:MODE',
        'acq_type': ':ACQ:TYPE',
        'acq_count': ':ACQ:AVER',
        'acq_memory_depth': ':ACQ:MDEP',
        'channel': ':WAV:SOUR',
        'mode': ':WAV:MODE',
        'format': ':WAV:FORM',
        'timebase': ':TIM:MAIN:SCAL',
        'offset': ':OFFS',
        'probe': ':PROB',
        'vert_scale': ':SCAL'
    }

    def __init__(self, name='None', settings=None):
        """

        Args:
            name (str): optional name of instance of class
            settings (list): list of other values to initialize class with

        """

        super(RigolOscilloscope, self).__init__(name, settings)

        # keep track of when the instrument was updated last to prevent sending requests to frequently
        self._last_update_time = time.time()

        rm = pyvisa.ResourceManager()

        # todo: JG 20170623 implement proper error handling when insturment is not connected.
        self.osci = rm.open_resource(self.settings['visa_resource'])
        self.osci.query_delay = 1000
        self._channel = self.settings['channel']
        baudrate = 9600

        self.osci.write_termination = '\n'
        self.osci.read_termination = '\n'

        self.osci.timeout = 10000
        # self.osci.inputbuffersize = 1000
        # self.osci.ByteOrder = 'littleEndian'
        self.osci.baud_rate = baudrate  #9600

        self.osci.write('*RST') #Places the oscilloscope in the factory default setup state.
        # self._wait_for_osci()
        # self._wait_for_osci()
        # self._wait_for_osci()
        # self._wait_for_osci()
        # print('waited')
        self.update(self.settings)
        # self.osci.write(':BUS1:RS232:BAUD ' + str(baudrate))
        # except:
        #     raise
    #
    # def __del__(self):
    #     self.osci.close()

    def _query(self, command):
        self.osci.write(command)
        return self.osci.read_raw().decode()[:-1]

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
          # Places the oscilloscope in the factory default setup state.

        super(RigolOscilloscope, self).update(settings)

        for key, value in settings.items():
            # print(key)
            if key not in ['visa_resource', 'connection_timeout']:
                command = self._COMMANDS[key]
                if key in ['offset', 'probe', 'vert_scale']:
                    command = ':' + self.settings['channel'] + command
                elif key in ['acq_type', 'acq_count', 'acq_memory_depth']:
                    self.osci.write(':SINGLE')
                self.osci.write(command + ' ' + str(value))

    def acq_time(self):
        """
        estimates the acquisition time
        Returns:

        """
        total_time = float(self.settings['timebase'])*10.0
        return total_time

    @staticmethod
    def _num_to_mdepth(mdepth):
        log = np.log10(mdepth)
        if log >= 6:
            return str(int(mdepth / 1e6)) + 'M'
        else:
            return str(int(mdepth / 1e3)) + 'k'

    def _set_acq(self):
        acq_params = ['acq_type', 'acq_count', 'acq_memory_depth']
        for param in acq_params:
            self.osci.write(self._COMMANDS[param] + ' ' + str(self.settings[param]))

    @property
    def _PROBES(self):
        return{'timebase_format':'time base format',
               'acq_type':'acquisition type',
               'acq_count':'acquisition count',
               'acq_memory_depth': 'number of points per trace',
               'mode': 'waveform mode',
               'format': 'waveform_format',
               'channel': 'waveform channel',
               'timebase': 'waveform timebase in seconds',
               'offset': 'waveform y offset',
               'probe': 'y intercept',
               'vert_scale': 'vertical scale per division in volts',
        }

    def read_probes(self, probe_name):

        # self._wait_for_osci()
        # assert(False)
        assert self._settings_initialized  # will cause read_probes to fail if settings (and thus also connection) not yet initialized
        if probe_name == 'visa_resource':
            return self.settings['visa_resource']
        if probe_name == 'connection_timeout':
            return self.settings['connection_timeout']

        assert probe_name in list(self._PROBES.keys())

        command = self._COMMANDS[probe_name] + '?'

        if probe_name in ['offset', 'probe', 'vert_scale']:
            channel = self.settings['channel']
            command = ':' + channel + command

        output = self._query(command)

        if probe_name == 'acq_memory_depth':
            return self._num_to_mdepth(float(output))
        elif probe_name in ['acq_count', 'probe']:
            return int(output)
        elif probe_name in ['timebase', 'offset', 'vert_scale']:
            return float(output)
        else:
            return output



    @property
    def is_connected(self):
        """
        Checks if the instrument is connected.
        Returns: True if connected, False otherwise.

        """
        identification = self.osci.query("*IDN?", delay = 0.1)
        return identification.strip("\n") == self._INSTRUMENT_IDENTIFIER

    def set_auto_trigger(self, timebase):
        self.osci.write(self._COMMANDS['timebase_format'] + ' MAIN')
        self.osci.write(self._COMMANDS['timebase'] + ' ' + str(timebase))
        self.osci.write(':TRIG:SWEep AUTO')
        self.osci.write(':RUN')

    def get_std_voltage(self):
        return float(self._query(':MEAS:ITEM? VAR, ' + self.settings['channel']))


    def get_timetrace(self):
        """
        reads a time trace from the intruments
        Returns:

        """
        half_wait = 10 * self.settings['timebase'] + 1
        self.osci.write(':TIM:MAIN:OFFS ' + str(5 * self.settings['timebase']))
        self.osci.write(':SINGLE') # start a single acquisition
        triggered = False
        while not triggered:
            time.sleep(half_wait)
            if self.osci.query('TRIG:STAT?', delay = 0.1) != 'RUN':
                triggered = True
        print('triggering')
        self.osci.write(':TFORce') # force trigger to make sure that the acquisition starts
        total_time = int(np.ceil(half_wait))
        for i in range(total_time):
            print(('waiting {:d}/{:d}'.format(i + 1, total_time)))
            time.sleep(1)

        # convert into dictionary
        # preamble = {k:float(v) for k, v in zip(['format', 'type', 'points', 'count', 'xincrement', 'xorigin', 'xreference', 'yincrement', 'yorigin', 'yreference'], preambleBlock.split(','))}
        yoff = int(self._query(':WAV:YOR?'))
        yzero = int(self._query(':WAV:YREF?'))
        ymult = float(self._query(':WAV:YINC?'))

        # self.osci.baud_rate = 18000000  #9600

        format = str(self.settings['format']).lower()
        # depending on the setting we get differnet data back
        if format == 'ascii':
            print('datatype ascii')
            print('JG: WARNING ascii NOT TESTED!')

            # send command to read data
            # raw_data = self.osci.query(':WAV:DATA?')
            # raw_data = self.osci.query_ascii_values(':WAV:DATA?', delay = 60)
            raw_data = self.osci.query(':WAV:DATA?', delay=120)
            # the first charactars are some metadata
            prefix = raw_data[:11]
            data = raw_data[11:].split(',')[:-1]
            # data = [float(d) for d in data]
            data = np.array(data, dtype=np.float32)
            # data = [(a - yzero) * ymult + yoff for a in data]


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
            self.osci.write(':WAV:DATA?')
            # # time.sleep(5)
            raw_data = self.osci.read_raw()
            data = [byte for byte in raw_data[11:]]
            data = [(a - yzero) * ymult + yoff for a in data][:-1]
        else:
            raise Exception('Unknown format')

        return data


    # def __del__(self):
    #     #COMMENT_ME
    #     # self._wait_for_osci()
    #     # self._set_mode('SpectrumAnalyzer')
    #     self.osci.close()
    #
    def _wait_for_osci(self):
        # #COMMENT_ME
        print('waiting')
        if self._last_update_time - time.time() < 1.0:
            print('sleeping')
            time.sleep(1)

        self._last_update_time = time.time()
        print('-- waited')
        # return


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

# class keysight_Oscilloscope(Oscilloscope):
#     """
#     This class provides a python implementation of the Keysight DSO1024A oscilloscope.
#     """
#
#     _INSTRUMENT_IDENTIFIER = 'Agilent Technologies,DSO1024A,CN56301388,00.04.06 SP06'
#     # String returned by spectrum analyzer upon querying it with '*IDN?'
#
#     _DEFAULT_SETTINGS = Parameter([
#         Parameter('visa_resource', 'USB0::0x0957::0x0588::CN56301388::0::INSTR', (str),
#                       'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
#         Parameter('timebase',[
#             Parameter('format', 'yt', ['yt', 'xy', 'roll'], 'time base format')
#         ]),
#         # Parameter('frequency_step', 10e6, float, 'frequency interval of spectrum analyzer frequency range'),
#         Parameter('acquisition',[
#                       Parameter('type','normal',['normal'], 'acquisition type'),
#                       Parameter('count',3, int, 'acquisition count')
#                   ]),
#         Parameter('waveform', [
#             Parameter('mode', 'raw', ['raw', 'normal'], 'waveform mode (not that normal allows only to acquire up to 600 datapoints! raw up tp 10240)'),
#             Parameter('points', 10240, int, 'waveform length, max is 10240 (mode raw) or 600 (mode normal)'),
#             Parameter('format', 'word', ['word', 'ascii', 'byte'], 'waveform format '),
#             Parameter('channel', 1, [1, 2, 3, 4], 'channel from which to read the data'),
#             Parameter('timebase', 1, [1,2,5, 10, 20, 50, 100, 200, 500], 'timebase: units per devision'),
#             Parameter('timebase_unit', 'ms', ['ns', 'us','ms', 's'], 'timebase: units per devision'),
#             Parameter('vert_scale', '5E-2', ['2E-3', '5E-3', '1E-2', '2E-2', '5E-2', '1E-1', '2E-1', '5E-1'], 'vert_scale volts per division (2 mV to 10 V)'),
#             #Parameter('offset', 0.0, float, 'vertical scale off set in V'),
#             Parameter('offset', '-1E-2', ['-5E-1', '-2E-1', '-1E-1', '-5E-2', '-2E-2', '-1E-2', '-5E-3', '-2E-3', '2E-3', '5E-3', '1E-2', '2E-2', '5E-2', '1E-1', '2E-1', '5E-1'],
#                       'offset on channel in volts'),
#             Parameter('probe', '1X', ['1X', '10X'], 'multipliation factor for y scale')
#         ]),
#         Parameter('trigger', [
#             Parameter('on', True, bool, 'trigger on or off'),
#             Parameter('channel', 2, [1, 2, 3, 4], 'channel from which to trigger'),
#             Parameter('mode', 'edge_pos', ['edge_pos', 'edge_neg'], 'trigger mode (psotive edge, negative edge)'),
#             Parameter('level', 0.0, [0.0, 0.1, 1], 'trigger level (I think in fractions of a devision between -6 and 6 could also be in volt)'),
#
#         ]),
#         Parameter('connection_timeout', 1000, int, 'the time to wait for a response from the oscilloscope with each query (units??)'),
#         ])
#
#
#
#     _PROBES = {}
#     # _PROBES = {'start_frequency': 'the lower bound of the frequency sweep',
#     #            'stop_frequency': 'the upper bound of the frequency sweep',
#     #            'trace': 'the frequency sweep of the inputted signal',
#     #            'tracking_generator': 'checks if the tracking generator is on',
#     #            'bandwidth': 'the curent bandwidth of the spectrum analyzer',
#     #            'output_power': 'the power of the tracking generator',
#     #            'mode': 'Spectrum Analyzer Mode or Tracking Generator Mode'}
#
#     def __init__(self, name='Oscilloscope', settings={}):
#         """
#
#         Args:
#             name (str): optional name of instance of class
#             settings (list): list of other values to initialize class with
#
#         """
#
#         super(Oscilloscope, self).__init__(name, settings)
#
#         # keep track of when the instrument was updated last to prevent sending requests to frequently
#         self._last_update_time = time.time()
#
#         rm = visa.ResourceManager()
#
#         # todo: JG 20170623 implement proper error handling when insturment is not connected.
#         self.osci = rm.open_resource(self.settings['visa_resource'])
#         self.osci.read_termination = '\n'
#         self.osci.timeout = self.settings['connection_timeout']
#         self.osci.inputbuffersize = 1000
#         self.osci.ByteOrder = 'littleEndian'
#
#         self.osci.write('*RST') #Places the oscilloscope in the factory default setup state.
#         self._wait_for_osci()
#         self.update(self.settings)
#         # except:
#         #     raise
#
#     def reset(self):
#         """
#         resets the oscilloscope to the default values
#         Returns:
#
#         """
#         self.osci.write('*RST') #Places the oscilloscope in the factory default setup state.
#
#     def update(self, settings):
#         """
#         updates the instrument parameters
#
#         Args:
#             settings: dictionary that contains the parameter indentifiers (keys) and the new parameters values (value)
#
#         """
#         super(Oscilloscope, self).update(settings)
#
#         if 'timebase' in settings:
#             # self._wait_for_osci()
#             # self._set_timebase(settings['timebase'])
#             self.osci.write(':TIMEBASE:MODE MAIN')
#             self.osci.write(':TIMEBASE:FORM ' + settings['timebase']['format'].upper())
#
#         if 'acquisition' in settings:
#             # self._wait_for_osci()
#             # self._set_acquisition(settings['acquisition'])
#             self.osci.write(':ACQUIRE:TYPE ' + settings['acquisition']['type'].upper())
#             self.osci.write(':ACQUIRE:COUNT ' + str(settings['acquisition']['count']))
#
#         if 'waveform' in settings:
#             channel = str(settings['waveform']['channel'])
#             waveform = settings['waveform']
#             # self._wait_for_osci()
#             # self._set_waveform(settings['waveform'])
#             self.osci.write(':WAV:POINTS:MODE ' + waveform['mode'].upper())
#             self.osci.write(':WAV:POINTS ' + str(waveform['points']))
#             self.osci.write(':WAV:FORMAT ' + waveform['format'].upper())
#             self.osci.write(':WAV:SOURCE CHAN' + channel)
#             self.osci.write(':TIM:SCAL ' + self.time_base_to_nr3(waveform['timebase'], waveform['timebase_unit']))
#             self.osci.write(':CHAN' + channel + ':PROB ' + waveform['probe']) #PROBE MUST BE SET BEFORE SCALE
#             self.osci.write(':CHAN' + channel + ':SCAL ' + waveform['vert_scale'])
#             self.osci.write(':CHAN' + channel + ':OFFS ' + waveform['offset'])#'{:0.02e}'.format(waveform['offset']))
#
#         if 'trigger' in settings:
#             channel = str(settings['trigger']['channel'])
#             self.osci.write(':TRIG:PULS:SOUR ' + channel)
#             if 'mode' in settings['trigger']:
#                 if settings['trigger'] is 'edge_pos':
#                     self.osci.write(':TRIG:EDGE:SLOP POS')
#                 if settings['trigger'] is 'edge_neg':
#                     self.osci.write(':TRIG:EDGE:SLOP NEG')
#
#             self.osci.write(':TRIG:EDGE:LEV ' + str(settings['trigger']['level']))
#
#
#     def acq_time(self):
#         """
#         estimates the acquisition time
#         Returns:
#
#         """
#         waveform = self.settings['waveform']
#         total_time = float(self.time_base_to_nr3(waveform['timebase'], waveform['timebase_unit']))*10
#         return total_time
#
#
#     def read_probes(self, probe_name):
#         print('NOT IMPLEMENTED')
#         pass
#         # self._wait_for_osci()
#         #
#         # if probe_name == 'start_frequency':
#         #     return self._get_start_frequency()
#         # elif probe_name == 'stop_frequency':
#         #     return self._get_stop_frequency()
#         # elif probe_name == 'trace':
#         #     return self._get_trace()
#         # elif probe_name == 'output_on':
#         #     return self._is_output_on()
#         # elif probe_name == 'bandwidth':
#         #     return self._get_bandwidth()
#         # elif probe_name == 'output_power':
#         #     return self._get_output_power()
#         # elif probe_name == 'mode':
#         #     return self._get_mode()
#         # else:
#         #     message = 'no probe with that name exists!'
#         #     raise AttributeError(message)
#
#     def is_connected(self):
#         """
#         Checks if the instrument is connected.
#         Returns: True if connected, False otherwise.
#
#         """
#         identification = self.osci.query('*IDN?')
#         return identification == self._INSTRUMENT_IDENTIFIER
#
#     def get_timetrace(self):
#         """
#         reads a time trace from the intruments
#         Returns:
#
#         """
#         self.osci.write(':SINGLE') # start a single acquisition
#         self.osci.write(':FORCetrig') # force trigger to make sure that the acquisition starts
#
#
#         # estimate acquisition time and wait until acquisition is finished
#         # JG: This seems to work, just wait enough time so that the oscilloscope can acquire the data
#         # output counter to terminal to show that the program didn't freeze
#         total_time = self.acq_time()
#         total_time = int(total_time )*10 # multiply by 10 because the wait loop time is 100ms
#         total_time = int(total_time * 1.2+1) # give another percent margin, this is empirical
#
#         for i in range(total_time):
#             print(('waiting {:d}/{:d}'.format(i, total_time)))
#             time.sleep(0.1)
#
#         operationComplete = bool(self.osci.query('*OPC?'))
#         print(('operationComplete', operationComplete))
#
#
#         # Get the preamble block
#         preambleBlock = self.osci.query(':WAV:PREAMBLE?')
#         # preable contains the curren settings
#         #   FORMAT        : int16 - 0 = WORD, 1 = BYTE, 2 = ASCII.
#         #   TYPE          : int16 - 0 = NORMAL, 1 = PEAK DETECT, 2 = AVERAGE
#         #   POINTS        : int32 - number of data points transferred.
#         #   COUNT         : int32 - 1 and is always 1.
#         #   XINCREMENT    : float64 - time difference between data points.
#         #   XORIGIN       : float64 - always the first data point in memory.
#         #   XREFERENCE    : int32 - specifies the data point associated with x-origin.
#         #   YINCREMENT    : float32 - voltage diff between data points.
#         #   YORIGIN       : float32 - value is the voltage at center screen.
#         #   YREFERENCE    : int32 - specifies the data point where y-origin occurs
#
#         # convert into dictionary
#         preamble = {k:float(v) for k, v in zip(['format', 'type', 'points', 'count', 'xincrement', 'xorigin', 'xreference', 'yincrement', 'yorigin', 'yreference'], preambleBlock.split(','))}
#
#
#
#         format = str(self.settings['waveform']['format']).lower()
#         # depending on the setting we get differnet data back
#         if format == 'ascii':
#             print('datatype ascii')
#             print('JG: WARNING ascii NOT TESTED!')
#
#             # send command to read data
#             # raw_data = self.osci.query(':WAV:DATA?')
#             raw_data = self.osci.query_ascii_values(':WAV:DATA?')
#             # the first charactars are some metadata
#             prefix = raw_data[:11]
#             data = raw_data[11:].split(',')
#             # data = [float(d) for d in data]
#             data = np.array(data)
#
#         elif format == 'word':
#             print('datatype word')
#             # try until we get the data
#             max_attempts = 100
#             count = 0
#             raw_data = []
#             while len(raw_data) == 0 and count<max_attempts:
#                 raw_data = self.osci.query_binary_values(':WAV:DATA?', datatype='H')
#                 print(('data empty retry, attempt {:02d}/{:02d}'.format(count, max_attempts)))
#                 count += 1
#                 time.sleep(0.1)
#             data = raw_data
#
#         elif format == 'byte':
#             print('datatype byte')
#             print('JG: WARNING BYTE NOT TESTED!')
#             raw_data = self.osci.query_binary_values(':WAV:DATA?')
#             data = raw_data
#         else:
#             print('WARNING UNKNOWN DATA FORMAT')
#
#         if len(data)>0:
#             dt = float(self.time_base_to_nr3(self.settings['waveform']['timebase'], self.settings['waveform']['timebase_unit']))*10/len(data)
#         else:
#             dt = 1
#         # add more meta data to the preamble
#         preamble['dt']= dt
#         preamble['vert_scale'] =float(self.settings['waveform']['vert_scale'])*float(self.settings['waveform']['probe'].split('X')[0]) #  vertical scale of osci
#         return data, preamble
#
#
#     def __del__(self):
#         #COMMENT_ME
#         # self._wait_for_osci()
#         # self._set_mode('SpectrumAnalyzer')
#         self.osci.close()
#
#     def _wait_for_osci(self):
#         #COMMENT_ME
#         print('waiting')
#         if self._last_update_time - time.time() < 1.0:
#             print('sleeping')
#             time.sleep(1)
#
#         self._last_update_time = time.time()
#         print('-- waited')
#
#
#     @staticmethod
#     def time_base_to_nr3(n, s):
#         """
#         converts the time base into the nr3 format
#         Args:
#             n: time base value (1, 2, 5)
#             s: time base (ns, us, ms, s)
#
#         Returns:
#             a string in nr3 format that represents the time for example, 1.0E-9, 2.0E-9, 5.0E-9, ... 1.0E+00, 2.0E+00, 5.0E+00)
#
#         """
#
#         if s == 'ns':
#             nr3 = 'E-9'
#         if s == 'us':
#             nr3 = 'E-6'
#         if s == 'ms':
#             nr3 = 'E-3'
#         if s == 's':
#             nr3 = 'E+00'
#
#         return str(n) + '.0' + nr3

if __name__ == '__main__':

        print('create oscilloscope instance:')
        oscil = rigol_Oscilloscope()
        print(oscil.is_connected())


        print('=============')
        # oscil.update({'waveform':{'CHAN1':{'vert_scale':0.5}}})

        # print((oscil.settings))

        # oscil.osci.write(':SINGLE')  # start a single acquisition
        # oscil.osci.write(':TFORce')

        data, preambleBlock = oscil.get_timetrace()

        dt = preambleBlock['dt']
        time = dt*np.arange(len(data))
        print(('data', data))
        plt.plot(time, data, '-x')

        plt.show()

