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

# class Oscilloscope(Instrument):

class rigol_Oscilloscope(Instrument):
    """
    This class provides a python implementation of the Keysight DSO1024A oscilloscope.
    """

    # String returned by spectrum analyzer upon querying it with '*IDN?'
    _INSTRUMENT_IDENTIFIER = 'RIGOL TECHNOLOGIES,MSO5104,MS5A241503457,00.01.03.00.03'

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_resource', 'USB0::0x1AB1::0x0515::MS5A241503457::INSTR', (str),
                      'pyVisa instrument identifier, to make a connection using the pyVisa package.'),
        Parameter('timebase',[
            Parameter('format', 'main', ['main', 'xy', 'roll'], 'time base format')
        ]),
        # Parameter('frequency_step', 10e6, float, 'frequency interval of spectrum analyzer frequency range'),
        Parameter('acquisition',[
                      Parameter('type','normal',['normal', 'averages', 'peak', 'hresolution'], 'acquisition type'),
                      Parameter('count',4.0, float, 'acquisition count')
                  ]),
        Parameter('waveform', [
            Parameter('mode', 'raw', ['normal', 'maximum', 'raw'], 'waveform mode (not that normal allows only to acquire up to 1000 datapoints!)'),
            Parameter('points', 1000, int, 'waveform length, max is 10240 (mode raw) or 1000 (mode normal)'),
            Parameter('format', 'byte', ['word', 'ascii', 'byte'], 'waveform format '),
            Parameter('channel', 'CHAN1', ['CHAN1','CHAN2','CHAN3','CHAN4'], 'channel from which to read the data'),
            Parameter('timebase', 0.000001, float, 'timebase: units per devision'),
            Parameter('CHAN1', [
                Parameter('offset', 0.0, float, 'offset on channel in volts. range is +/- 1V, +/1 30V, +/- 100V depending on vertical scale'),
                Parameter('probe', 1.0, float,
                          'multipliation factor for y scale - can only be 1s, 2s, 5s from 1e-4 to 1e4'),
                Parameter('vert_scale', 100.E-3, float, 'volts per division (500 uV to 10 V)')]),
            Parameter('CHAN2', [
                Parameter('offset', 0.0, float,
                          'offset on channel in volts. range is +/- 1V, +/1 30V, +/- 100V depending on vertical scale'),
                Parameter('probe', 1.0, float,
                          'multipliation factor for y scale - can only be 1s, 2s, 5s from 1e-4 to 1e4'),
                Parameter('vert_scale', 100.E-3, float, 'volts per division (500 uV to 10 V)')]),
            Parameter('CHAN3', [
                Parameter('offset', 0.0, float,
                          'offset on channel in volts. range is +/- 1V, +/1 30V, +/- 100V depending on vertical scale'),
                Parameter('probe', 1.0, float,
                          'multipliation factor for y scale - can only be 1s, 2s, 5s from 1e-4 to 1e4'),
                Parameter('vert_scale', 100.E-3, float, 'volts per division (500 uV to 10 V)')]),
            Parameter('CHAN4', [
                Parameter('offset', 0.0, float,
                          'offset on channel in volts. range is +/- 1V, +/1 30V, +/- 100V depending on vertical scale'),
                Parameter('probe', 1.0, float,
                          'multipliation factor for y scale - can only be 1s, 2s, 5s from 1e-4 to 1e4'),
                Parameter('vert_scale', 100.E-3, float, 'volts per division (500 uV to 10 V)')]),
        ]),
        Parameter('trigger', [
            # Parameter('on', True, bool, 'trigger on or off'),
            Parameter('channel','CHAN1', ['CHAN1','CHAN2','CHAN3','CHAN4'], 'channel from which to trigger'),
            Parameter('mode', 'POS', ['POS', 'NEG', 'RFAL'], 'trigger mode (psotive edge, negative edge)'),
            Parameter('level', 0.0, float, 'trigger level, range is (-5 x vertical scale - offset) to (5 X vertical scale - offset)'),

        ]),
        Parameter('connection_timeout', 9., float, 'the time to wait for a response from the oscilloscope with each query in seconds (16 ns to 10s'),
        ])

    def __init__(self, name='Oscilloscope', settings={}):
        """

        Args:
            name (str): optional name of instance of class
            settings (list): list of other values to initialize class with

        """

        super(rigol_Oscilloscope, self).__init__(name, settings)

        # keep track of when the instrument was updated last to prevent sending requests to frequently
        self._last_update_time = time.time()

        rm = visa.ResourceManager()

        print(rm.list_resources())

        # todo: JG 20170623 implement proper error handling when insturment is not connected.
        self.osci = rm.open_resource(self.settings['visa_resource'])

        self.osci.read_termination = '\n'
        self.osci.write_termination = '\n'
        self.osci.timeout = self.settings['connection_timeout']
        self.osci.inputbuffersize = 1000
        self.osci.ByteOrder = 'littleEndian'
        self.osci.baud_rate = 9600
        #
        self.osci.write('*RST') #Places the oscilloscope in the factory default setup state.
        self._wait_for_osci()
        self._wait_for_osci()
        self._wait_for_osci()
        self._wait_for_osci()
        print('waited')
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
        super(rigol_Oscilloscope, self).update(settings)

        for key, value in settings.items():
            # print(key)
            if key == 'visa_resource' or key == 'connection_timeout':
                continue
            for subkey, subvalue in value.items():
                # print(subkey)
                # print(self.is_connected())
                subvalue = self._value_to_internal(subvalue)
                if key == 'timebase':
                    subkey = self._timebase_param_to_internal(subkey)
                elif key == 'acquisition':
                    subkey = self._acq_param_to_internal(subkey)
                elif key == 'waveform':
                    if subkey in ['CHAN1','CHAN2','CHAN3','CHAN4']:
                        subkey = self._wav_param_to_internal(subkey, optional_dict = self.settings['waveform'][subkey])
                        continue
                    else:
                        subkey = self._wav_param_to_internal(subkey)
                elif key == 'trigger':
                    subkey = self._trig_param_to_internal(subkey)

                if self.is_connected():
                    print(subkey + ' ' + str(subvalue))
                    self.osci.write(subkey + ' ' + str(subvalue))
                    self._wait_for_osci()
                    self._wait_for_osci()


    def _value_to_internal(self, value):
        return str(value)

    def _timebase_param_to_internal(self, param):
        if param == 'format':
            return ':TIME:MODE'
    def _acq_param_to_internal(self, param):
        if param == 'type':
            return ':ACQ:TYPE'
        if param == 'count':
            return ':ACQ:AVER'
        else:
            raise KeyError

    def _wav_param_to_internal(self, param, optional_dict = None):
        if param in ['CHAN1', 'CHAN2', 'CHAN3', 'CHAN4']:
            print(optional_dict)
            for key, val in optional_dict.items():
                if key == 'offset':
                    if self.is_connected():
                        print(':'+param + ':OFFS ' + str(val))
                        self.osci.write(':'+param + ':OFFS ' + str(val))
                        self._wait_for_osci()
                    # return ':' + param + ':OFFS'
                elif key == 'probe':
                    if self.is_connected():
                        print(':'+param + ':PROB ' + str(val))
                        self.osci.write(':'+param + ':PROB ' + str(val))
                        self._wait_for_osci()
                    # return ':' + param + ':PROB'
                elif key == 'vert_scale':
                    if self.is_connected():
                        print(':'+param + ':SCAL ' + str(val))
                        self.osci.write(':'+param + ':SCAL ' + str(val))
                        self._wait_for_osci()
                    # return ':' + param + ':SCAL'
                else:
                    raise KeyError
        elif param == 'mode':
            return ':WAV:MODE'
        elif param == 'points':
            return ':WAV:POINT'
        elif param == 'format':
            return ':WAV:FORM'
        elif param == 'channel':
            return ':WAV:SOUR'
        elif param == 'timebase':
            return ':TIME:SCAL'
        else:
            raise KeyError

    def _trig_param_to_internal(self, param):
        if param == 'channel':
            return ':TRIG:PULS:SOUR'
        if param == 'mode':
            return ':TRIG:EDGE:SLOP'
        if param == 'level':
            return ':TRIG:EDGE:LEV'
        else:
            raise KeyError


        # if 'timebase' in settings:
        #     timebase = settings['timebase']
        #     # self._wait_for_osci()
        #     # self._set_timebase(settings['timebase'])
        #     # self.osci.write(':TIM:MODE MAIN')
        #     if
        #     self.osci.write(':TIM:MODE ' + settings['timebase']['format'])
        #
        # if 'acquisition' in settings:
        #     # self._wait_for_osci()
        #     # self._set_acquisition(settings['acquisition'])
        #     self.osci.write(':ACQ:TYPE ' + settings['acquisition']['type'])
        #     self.osci.write(':ACQ:AVER ' + str(settings['acquisition']['count']))
        #
        # if 'waveform' in settings:
        #     channel = str(settings['waveform']['channel'])
        #     waveform = settings['waveform']
        #     # self._wait_for_osci()
        #     # self._set_waveform(settings['waveform'])
        #     self.osci.write(':WAV:MODE ' + waveform['mode'])
        #     self.osci.write(':WAV:POINT ' + str(waveform['points']))
        #     self.osci.write(':WAV:FORM ' + waveform['format'])
        #     self.osci.write(':WAV:SOUR CHAN' + str(channel))
        #     self.osci.write(':TIM:SCAL ' + self.time_base_to_nr3(waveform['timebase'], waveform['timebase_unit']))
        #     self.osci.write(':CHAN' + str(channel) + ':PROB ' + str(waveform['probe'])) #PROBE MUST BE SET BEFORE SCALE
        #     self.osci.write(':CHAN' + str(channel) + ':SCAL ' + str(waveform['vert_scale']))
        #     self.osci.write(':CHAN' + str(channel) + ':OFFS ' + str(waveform['offset']))#'{:0.02e}'.format(waveform['offset']))
        #
        # if 'trigger' in settings:
        #     channel = str(settings['trigger']['channel'])
        #     self.osci.write(':TRIG:PULS:SOUR CHAN' + str(channel))
        #     if 'mode' in settings['trigger']:
        #         if settings['trigger'] == 'POS':
        #             self.osci.write(':TRIG:EDGE:SLOP POS')
        #         if settings['trigger'] == 'NEG':
        #             self.osci.write(':TRIG:EDGE:SLOP NEG')
        #         if settings['trigger'] == 'RFAL':
        #             self.osci.write(':TRIG:EDGE:SLOP RFAL')
        #
        #     self.osci.write(':TRIG:EDGE:LEV ' + str(settings['trigger']['level']))


    def acq_time(self):
        """
        estimates the acquisition time
        Returns:

        """
        waveform = self.settings['waveform']
        total_time = float(waveform['timebase'])
        return total_time

    def _PROBES(self):
        return{'timebase_format':'time base format',
               'acquisition_type':'acquisition type',
               'acquisition_count':'acquisition count',
               'waveform_mode':'waveform mode',
               'num_pts':'number of points',
               'waveform_format':'waveform_format',
               'waveform_channel':'waveform channel',
               'waveform_timebase':'waveform timebase in seconds',
               'trigger_channel':'channel from which to trigger',
               'trigger_mode':'trigger mode',
               'trigger_level':'trigger level'
        }

    def read_probes(self, probe_name):

        self._wait_for_osci()

        assert (self._settings_initialized)  # will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        if prob_name == 'timebase_format':
            return float(self.osci.query(':TIM:MODE?'))
        elif prob_name == 'acquisition_type':
            return self.osci.query(':ACQ:TYPE?')
        elif prob_name == 'acquisition_count':
            return int(self.osci.query(':ACQ:AVER?'))
        elif prob_name == 'waveform_mode':
            return self.osci.query(':WAV:MODE?')
        elif prob_name == 'num_pts':
            return int(self.osci.query(':WAV:POIN?'))
        elif prob_name == 'waveform_format':
            return self.osci.query(':WAV:FORM?')
        elif prob_name == 'waveform_channel':
            return self.osci.query(':WAV:SOUR?')
        elif prob_name == 'waveform_timebase':
            return float(self.osci.query(':TIME:MAIN:SCAL?'))
        elif prob_name == 'trigger_channel':
            return self.osci.query(':TRIG:PULS:SOUR?')
        elif prob_name == 'trigger_mode':
            return self.osci.query(':TRIG:EDGE:SLOP?')
        elif prob_name == 'trigger_level':
            return self.osci.query(':TRIG:EDGE:LEV?')
        else:
            raise KeyError


    def is_connected(self):
        """
        Checks if the instrument is connected.
        Returns: True if connected, False otherwise.

        """
        identification = self.osci.query('*IDN?')
        # self.osci.write('*IDN?')
        # time.sleep(0.1)
        # print(self.osci.read())
        return identification.strip("\n") == self._INSTRUMENT_IDENTIFIER

    def get_timetrace(self):
        """
        reads a time trace from the intruments
        Returns:

        """
        self.osci.write(':SINGLE') # start a single acquisition
        self._wait_for_osci()
        self.osci.write(':TFORce') # force trigger to make sure that the acquisition starts
        self._wait_for_osci()
        self.is_connected()

        # estimate acquisition time and wait until acquisition is finished
        # JG: This seems to work, just wait enough time so that the oscilloscope can acquire the data
        # output counter to terminal to show that the program didn't freeze
        total_time = self.acq_time()
        total_time = int(total_time )*10 # multiply by 10 because the wait loop time is 100ms
        total_time = int(total_time * 1.2+1) # give another percent margin, this is empirical

        for i in range(total_time):
            print(('waiting {:d}/{:d}'.format(i, total_time)))
            time.sleep(0.1)
        if self.is_connected():
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
        self._wait_for_osci()
        print(preamble)

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
            print(self.is_connected())
            raw_data = self.osci.query_binary_values(':WAV:DATA?', delay = 3)
            # raw_data = self.osci.query(':WAV:DATA?')
            data = raw_data
        else:
            print('WARNING UNKNOWN DATA FORMAT')

        if len(data)>0:
            dt = float((self.settings['waveform']['timebase'])*10/len(data))
        else:
            dt = 1
        # add more meta data to the preamble
        preamble['dt']= dt
        channel = self.settings['waveform']['channel']
        # preamble['vert_scale'] =float(self.settings['waveform'][channel]['vert_scale'])*float(self.settings['waveform'][channel]['probe'].split('X')[0]) #  vertical scale of osci
        preamble['vert_scale'] = float(self.settings['waveform'][channel]['vert_scale']) * float(self.settings['waveform'][channel]['probe'])  # vertical scale of osci
        return data, preamble


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

        # oscil.write(':SINGLE')  # start a single acquisition
        # oscil.write(':TFORce')

        data, preambleBlock = oscil.get_timetrace()

        dt = preambleBlock['dt']
        time = dt*np.arange(len(data))
        print(('data', data))
        plt.plot(time, data, '-x')

        plt.show()
