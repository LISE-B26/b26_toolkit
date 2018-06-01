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

#TODO:
# tasklist setter


import ctypes
import os
import numpy
import warnings
from pylabcontrol.core.read_write_functions import get_config_value

from pylabcontrol.core import Instrument, Parameter
import numpy as np

##############################
# Setup some typedefs and constants
# to correspond with values in
# C:\Program Files\National Instruments\NI-DAQ\DAQmx ANSI C Dev\include\NIDAQmx.h
# or http://digital.ni.com/public.nsf/ad0f282819902a1986256f79005462b1/b77ebfb849f162cd86256f150048dbb1/$FILE/NIDAQmx.h
# the typedefs
int32 = ctypes.c_long
int64 = ctypes.c_longlong
uInt8 = ctypes.c_uint8
uInt32 = ctypes.c_ulong
uInt64 = ctypes.c_ulonglong
float64 = ctypes.c_double
bool32 = ctypes.c_bool
TaskHandle = uInt64
# Analog constants
DAQmx_Val_Cfg_Default = int32(-1)
DAQmx_Val_Volts = 10348
DAQmx_Val_Rising = 10280
DAQmx_Val_FiniteSamps = 10178
DAQmx_Val_ContSamps = 10123
DAQmx_Val_GroupByChannel = 0

# DI constants
DAQmx_Val_CountUp = 10128
DAQmx_Val_Hz = 10373  # Hz
DAQmx_Val_Low = 10214  # Low
DAQmx_Val_Seconds = 10364
DAQmx_Val_Ticks = 10304  # specifies units as timebase ticks

DAQmx_Val_ChanPerLine = 0  # One Channel For Each Line
DAQmx_Val_ChanForAllLines = 1  # One Channel For All Lines


# =============== NI DAQ 6259======= =======================
# ==========================================================

class DAQ(Instrument):
    """
    Class containing all functions used to interact with the NI DAQ, mostly
    acting as a wrapper around C-level dlls provided by NI. Tested on an
    NI DAQ 6259, but should be compatable with most daqmx devices. Supports
    analog output (ao), analog input (ai), and digital input (di) channels.
    Also supports gated digital input, using one PFI channel as a counter
    and a second as a clock.

    In general, the order of calls to use one of these channels is:
    setup
    run
    (read or write if required)
    stop

    In general, 'setup' sets up the buffer, either filling it with the values to
    output or telling it to ready for input, and locks the task so the correct clock.
    'Run' starts the input to or output from the buffer.
    'Read' sends data from the buffer to the computer (if applicable).
    'Stop' ends the task and cleans up.
    """

    try:
        dll_path = get_config_value('NIDAQ_DLL_PATH',
                                    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
        if os.name == 'nt':
            # checks for windows. If not on windows, check for your OS and add
            # the path to the DLL on your machine
            nidaq = ctypes.WinDLL(dll_path)  # load the DLL
            dll_detected = True
        else:
            warnings.warn("NI DAQmx DLL not found. If it should be present, check the path:")
            print(dll_path)
            dll_detected = False
    except WindowsError:
        # make a fake DAQOut instrument
        dll_detected = False
    except:
        raise

    tasklist = {}
    tasknum = 0

    # currently includes four analog outputs, five analog inputs, and one digital counter input. Add
    # more as needed and your device allows
    _DEFAULT_SETTINGS = Parameter([
        Parameter('device', 'Dev1', ['Dev1', 'cDAQ9184-1BA7633Mod3', 'cDAQ9184-1BA7633Mod4', 'cDAQ9184-1BA7633Mod1', 'cDAQ1Mod1'], 'Name of DAQ device'),
        Parameter('override_buffer_size', -1, int, 'Buffer size for manual override (unused if -1)'),
        Parameter('ao_read_offset', .005, float, 'Empirically determined offset for reading ao voltages internally'),
        Parameter('analog_output',
                  [
                      Parameter('ao0',
                                [
                                    Parameter('channel', 0, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao1',
                                [
                                    Parameter('channel', 1, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao2',
                                [
                                    Parameter('channel', 2, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao3',
                                [
                                    Parameter('channel', 3, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                )
                  ]
                  ),
        Parameter('analog_input',
                  [
                      Parameter('ai0',
                                [
                                    Parameter('channel', 0, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai1',
                                [
                                    Parameter('channel', 1, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai2',
                                [
                                    Parameter('channel', 2, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai3',
                                [
                                    Parameter('channel', 3, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai4',
                                [
                                    Parameter('channel', 4, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage (V)')
                                ]
                                )
                  ]
                  ),
        Parameter('digital_input',
                  [
                      Parameter('ctr0',
                                [
                                    Parameter('input_channel', 0, list(range(0, 32)), 'channel for counter signal input'),
                                    Parameter('counter_PFI_channel', 8, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('gate_PFI_channel', 14, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('clock_PFI_channel', 13, list(range(0, 32)), 'PFI for clock channel output'),
                                    Parameter('clock_counter_channel', 1, [0, 1], 'channel for clock output'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)')
                                ]
                                ),
                      Parameter('ctr1',
                                [
                                    Parameter('input_channel', 1, list(range(0, 32)), 'channel for counter signal input'),
                                    Parameter('counter_PFI_channel', 3, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('gate_PFI_channel', 14, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('clock_PFI_channel', 12, list(range(0, 32)), 'PFI for clock channel output'),
                                    Parameter('clock_counter_channel', 0, [0, 1], 'channel for clock output'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)')
                                ]
                                )
                  ]
                  ),
        Parameter('digital_output',
                  [
                      Parameter('do0',
                                [
                                    Parameter('channel', 0, list(range(0, 16)), 'channel')
                                ]
                                ),
                      Parameter('do8',
                                [
                                    Parameter('channel', 8, list(range(0, 16)), 'channel')
                                ]
                                )
                  ]
                  )
    ])

    def __init__(self, name=None, settings=None):
        if self.dll_detected:
            # buf_size = 10
            # data = ctypes.create_string_buffer('\000' * buf_size)
            # try:
            #     #Calls arbitrary function to check connection
            #     self.CHK(self.nidaq.DAQmxGetDevProductType(device, ctypes.byref(data), buf_size))
            #     self.hardware_detected = True
            # except RuntimeError:
            #     self.hardware_detected = False
            super(DAQ, self).__init__(name, settings)

    def update(self, settings):
        """
        Updates daq settings for each channel in the software instrument.
        Unlike most instruments, all of the settings are sent to the DAQ on instantiation of
        a task, such as an input or output. Thus, changing the settings only updates the internal
        daq construct in the program and makes no hardware changes.
        Args:
            settings: a settings dictionary in the standard form
        """
        super(DAQ, self).update(settings)
        for key, value in settings.items():
            if key == 'device':
                if not (self.is_connected):
                    raise EnvironmentError('Device invalid, cannot connect to DAQ')

    def _add_to_tasklist(self, name, task):
        matching = [x for x in self.tasklist if name in x]
        if not matching:
            task_name = name + '000'
        else:
            last_task = sorted(matching)[-1]
            task_name = name + '{0:03d}'.format(int(last_task[-3:])+1)
        self.tasklist.update({task_name: task})
        return task_name

    @property
    def _PROBES(self):
        return None

    def read_probes(self, key):
        pass

    @property
    def is_connected(self):
        """
        Makes a non-state-changing call (a get id call) to check connection to a daq
        Returns: True if daq is connected, false if it is not
        """
        buf_size = 10
        data = ctypes.create_string_buffer(('\000' * buf_size).encode('ascii'))
        try:
            # Calls arbitrary function to check connection
            self._check_error(self.nidaq.DAQmxGetDevProductType(self.settings['device'].encode('ascii'), ctypes.byref(data), buf_size))
            return True
        except RuntimeError:
            return False

    def setup_counter(self, channel, sample_num, continuous_acquisition=False):
        """
        Initializes a hardware-timed digital counter, bound to a hardware clock
        Args:
            channel: digital channel to initialize for read in
            sample_num: number of samples to read in for finite operation, or number of samples between
                       reads for continuous operation (to set buffer size)
            continuous_acquisition: run in continuous acquisition mode (ex for a continuous counter) or
                                    finite acquisition mode (ex for a scan, where the number of samples needed
                                    is known a priori)

        Returns: source of clock that this method sets up, which can be given to another function to synch that
        input or output to the same clock

        """

        # Note that for this counter, we have two tasks. The normal 'task_handle' corresponds to the clock, and this
        # is the task which is started when run is called. The second 'task_handle_ctr' corresponds to the counter,
        # and this waits for the clock and will be started simultaneously.
        task = {
            'task_handle': None,
            'task_handle_ctr': None,
            'counter_out_PFI_str': None,
            'sample_num': None,
            'sample_rate': None,
            'num_samples_per_channel': None,
            'timeout': None
        }

        task_name = self._add_to_tasklist('ctr', task)

        if 'digital_input' not in list(self.settings.keys()):
            raise ValueError('This DAQ does not support digital input')
        if not channel in list(self.settings['digital_input'].keys()):
            raise KeyError('This is not a valid digital input channel')
        channel_settings = self.settings['digital_input'][channel]
        self.running = True
        task['sample_num'] = sample_num
        task['sample_rate'] = float(channel_settings['sample_rate'])
        if not continuous_acquisition:
            task['num_samples_per_channel'] = task['sample_num']
        else:
            task['num_samples_per_channel'] = -1
        task['timeout'] = float64(5 * (1 / task['sample_rate']) * task['sample_num'])
        input_channel_str = (self.settings['device'] + '/' + channel).encode('ascii')
        task['counter_out_PFI_str'] = ('/' + self.settings['device'] + '/PFI' + str(
            channel_settings['clock_PFI_channel'])).encode('ascii')  # initial / required only here, see NIDAQ documentation
        counter_out_str = (self.settings['device'] + '/ctr' + str(channel_settings['clock_counter_channel'])).encode('utf-8')
        task['task_handle_ctr'] = TaskHandle(0)
        task['task_handle'] = TaskHandle(1)

        # set up clock
        self._dig_pulse_train_cont(task, .5, counter_out_str)
        # set up counter using clock as reference
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle_ctr'])))
        self._check_error(self.nidaq.DAQmxCreateCICountEdgesChan(task['task_handle_ctr'],
                                                                 input_channel_str, "", DAQmx_Val_Rising, 0,
                                                                 DAQmx_Val_CountUp))
        # PFI13 is standard output channel for ctr1 channel used for clock and
        # is internally looped back to ctr1 input to be read
        if not continuous_acquisition:
            self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle_ctr'], task['counter_out_PFI_str'],
                                                               float64(task['sample_rate']), DAQmx_Val_Rising,
                                                               DAQmx_Val_FiniteSamps, uInt64(task['sample_num'])))
        else:
            self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle_ctr'], task['counter_out_PFI_str'],
                                                               float64(task['sample_rate']), DAQmx_Val_Rising,
                                                               DAQmx_Val_ContSamps, uInt64(task['sample_num'])))
        # if (self.settings['override_buffer_size'] > 0):
        #     self._check_error(self.nidaq.DAQmxCfgInputBuffer(self.DI_taskHandleCtr, uInt64(self.settings['override_buffer_size'])))
        # self._check_error(self.nidaq.DAQmxCfgInputBuffer(self.DI_taskHandleCtr, uInt64(sampleNum)))

        self._check_error(self.nidaq.DAQmxStartTask(task['task_handle_ctr']))

        return task_name

    def _dig_pulse_train_cont(self, task, DutyCycle, counter_out_str):
        """
        Initializes a digital pulse train to act as a reference clock
        Args:
            Freq: frequency of reference clock
            DutyCycle: percentage of cycle that clock should be high voltage (usually .5)
            Samps: number of samples to generate

        Returns:

        """
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle'])))
        self._check_error(self.nidaq.DAQmxCreateCOPulseChanFreq(task['task_handle'],
                                                                counter_out_str, '', DAQmx_Val_Hz, DAQmx_Val_Low,
                                                                float64(0.0),
                                                                float64(task['sample_rate']), float64(DutyCycle)))
        self._check_error(self.nidaq.DAQmxCfgImplicitTiming(task['task_handle'],
                                                            DAQmx_Val_ContSamps, uInt64(task['sample_num'])))

    def setup_gated_counter(self, channel, num_samples):
        """
        Initializes a gated digital input task. The gate acts as a clock for the counter, so if one has a fast ttl source
        this allows one to read the counter for a shorter time than would be allowed by the daq's internal clock.
        Args:
            channel: channel to use for counter input
            num_samples: number of samples to read on counter
        """
        if 'digital_input' not in list(self.settings.keys()):
            raise ValueError('This DAQ does not support digital input')
        if not channel in list(self.settings['digital_input'].keys()):
            raise KeyError('This is not a valid digital input channel')
        channel_settings = self.settings['digital_input'][channel]

        task = {
            'task_handle': None,
            'sample_num': None,
            'sample_rate': None,
            'num_samples_per_channel': None,
            'timeout': None
        }

        task_name = self._add_to_tasklist('gatedctr', task)

        input_channel_str_gated = (self.settings['device'] + '/' + channel).encode('ascii')
        counter_out_PFI_str_gated = ('/' + self.settings['device'] + '/PFI' + str(
            channel_settings['counter_PFI_channel'])).encode('ascii')  # initial / required only here, see NIDAQ documentation
        gate_PFI_str = ('/' + self.settings['device'] + '/PFI' + str(
            channel_settings['gate_PFI_channel'])).encode('ascii')  # initial / required only here, see NIDAQ documentation

        #set both to same value, no option for continuous counting (num_samples_per_channel == -1) with gated counter
        task['sample_num'] = num_samples
        task['num_samples_per_channel'] = num_samples

        task['task_handle'] = TaskHandle(0)

        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle'])))

        MIN_TICKS = 0;
        MAX_TICKS = 100000;


        # setup counter to measure pulse widths
        self._check_error(
            self.nidaq.DAQmxCreateCIPulseWidthChan(task['task_handle'], input_channel_str_gated, '', MIN_TICKS,
                                                   MAX_TICKS, DAQmx_Val_Ticks, DAQmx_Val_Rising, ''))

        # specify number of samples to acquire
        self._check_error(self.nidaq.DAQmxCfgImplicitTiming(task['task_handle'],
                                                            DAQmx_Val_FiniteSamps, uInt64(task['sample_num'])))

        # set the terminal for the counter timebase source to the APD source
        # in B26, this is the ctr0 source PFI8, but this will vary from daq to daq
        self._check_error(self.nidaq.DAQmxSetCICtrTimebaseSrc(task['task_handle'], input_channel_str_gated,
                                                              counter_out_PFI_str_gated))

        # set the terminal for the gate to the pulseblaster source
        # in B26, due to crosstalk issues when we use the default PFI9 which is adjacent to the ctr0 source, we set this
        # to the non-default value PFI14
        self._check_error(self.nidaq.DAQmxSetCIPulseWidthTerm(task['task_handle'], input_channel_str_gated,
                                                                  gate_PFI_str))

        # turn on duplicate count prevention (allows 0 counts to be a valid count for clock ticks during a gate, even
        # though the timebase never went high and thus nothing would normally progress, by also referencing to the internal
        # clock at max frequency, see http://zone.ni.com/reference/en-XX/help/370466AC-01/mxdevconsid/dupcountprevention/
        # for more details)
        self._check_error(
            self.nidaq.DAQmxSetCIDupCountPrevent(task['task_handle'], input_channel_str_gated, bool32(True)))

        return task_name

    # read sampleNum previously generated values from a buffer, and return the
    # corresponding 1D array of ctypes.c_double values
    def read_counter(self, task_name):
        """
        read sampleNum previously generated values from a buffer, and return the
        corresponding 1D array of ctypes.c_double values
        Returns: 1d array of ctypes.c_double values with the requested counts. Counts as given by the daq are a running
            total, that is if you get 5 counts/s, the returned array will be [5,10,15,20...]

        """
        task = self.tasklist[task_name]

        # difference between gated and non gated counter: the non-gated has a separate clock, while the gated one doesn't
        # For the gated case the task it self is also the clock
        # so if there is not extra handle for the clock we use the task_handle
        if 'task_handle_ctr' in task:
            task_handle_ctr = task['task_handle_ctr']
        else:
            task_handle_ctr = task['task_handle']

        # initialize array and integer to pass as pointers
        data = (float64 * task['sample_num'])()
        samplesPerChanRead = int32()


        self._check_error(self.nidaq.DAQmxReadCounterF64(task_handle_ctr,
                                                         int32(task['num_samples_per_channel']), float64(-1),
                                                         ctypes.byref(data),
                                                         uInt32(task['sample_num']),
                                                         ctypes.byref(samplesPerChanRead),
                                                         None))

        return data, samplesPerChanRead

    def setup_AO(self, channels, waveform, clk_source=""):
        """
        Initializes a arbitrary number of analog output channels to output an arbitrary waveform
        Args:
            channels: List of channels to output on
            waveform: 2d array of voltages to output, with each column giving the output values at a given time
                (the timing given by the sample rate of the channel) with the channels going from top to bottom in
                the column in the order given in channels
            clk_source: the PFI channel of the hardware clock to lock the output to, or "" to use the default
                internal clock
        """
        if 'analog_output' not in list(self.settings.keys()):
            raise ValueError('This DAQ does not support analog output')
        for c in channels:
            if not c in list(self.settings['analog_output'].keys()):
                raise KeyError('This is not a valid analog output channel')

        task = {
            'task_handle': None,
            'sample_num': None,
            'sample_rate': None,
            'num_samples_per_channel': None,
            'timeout': None
        }

        task_name = self._add_to_tasklist('ao', task)

        task['sample_rate'] = float(
            self.settings['analog_output'][channels[0]]['sample_rate'])  # float prevents truncation in division

        for c in channels:
            if not self.settings['analog_output'][c]['sample_rate'] == task['sample_rate']:
                raise ValueError('All sample rates must be the same')
        channel_list = ''.encode('ascii')
        for c in channels:
            channel_list += (self.settings['device'] + '/' + c + ',').encode('ascii')
        channel_list = channel_list[:-1]
        self.running = True
        # special case 1D waveform since length(waveform[0]) is undefined
        if (len(numpy.shape(waveform)) == 2):
            numChannels = len(waveform)
            task['sample_num'] = len(waveform[0])
        else:
            task['sample_num'] = len(waveform)
            numChannels = 1
        task['task_handle'] = TaskHandle(0)
        # special case 1D waveform since length(waveform[0]) is undefined
        # converts python array to ctypes array
        if len(numpy.shape(waveform)) == 2:
            data = numpy.zeros((numChannels, task['sample_num']),dtype=numpy.float64)
            for i in range(numChannels):
                for j in range(task['sample_num']):
                    data[i, j] = waveform[i, j]
        else:
            data = numpy.zeros((task['sample_num']), dtype=numpy.float64)
            for i in range(task['sample_num']):
                data[i] = waveform[i]

        if not (clk_source == ""):
            clk_source = self.tasklist[clk_source]['counter_out_PFI_str']

        self._check_error(self.nidaq.DAQmxCreateTask("",
                                                     ctypes.byref(task['task_handle'])))
        self._check_error(self.nidaq.DAQmxCreateAOVoltageChan(task['task_handle'],
                                                              channel_list,
                                                              "",
                                                              float64(-10.0),
                                                              float64(10.0),
                                                              DAQmx_Val_Volts,
                                                              None))
        self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle'],
                                                           clk_source,
                                                           float64(task['sample_rate']),
                                                           DAQmx_Val_Rising,
                                                           DAQmx_Val_FiniteSamps,
                                                           uInt64(task['sample_num'])))
        self._check_error(self.nidaq.DAQmxWriteAnalogF64(task['task_handle'],
                                                         int32(task['sample_num']),
                                                         0,
                                                         float64(-1),
                                                         DAQmx_Val_GroupByChannel,
                                                         data.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
                                                         None,
                                                         None))


        return task_name

    def setup_AI(self, channel, num_samples_to_acquire, continuous = False):
        """
        Initializes an input channel to read on
        Args:
            channel: Channel to read input
            num_samples_to_acquire: number of samples to acquire on that channel
        """

        task = {
            'task_handle': None,
            'sample_num': None,
            'sample_rate': None,
            'num_samples_per_channel': None,
            'timeout': None
        }

        task_name = self._add_to_tasklist('ai', task)

        channel_list = ''
        channel_list += self.settings['device'] + '/' + channel + ','

        if 'analog_input' not in list(self.settings.keys()):
            raise ValueError('This DAQ does not support analog input')
        task['task_handle'] = TaskHandle(0)
        task['sample_num'] = num_samples_to_acquire
        data = numpy.zeros((task['sample_num'],), dtype=numpy.float64)
        # now, on with the program
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle'])))
        self._check_error(self.nidaq.DAQmxCreateAIVoltageChan(task['task_handle'], channel_list, '',
                                                              DAQmx_Val_Cfg_Default,
                                                              float64(-10.0), float64(10.0),
                                                              DAQmx_Val_Volts, None))
        if not continuous:
            self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle'], "", float64(
                                                           self.settings['analog_input'][channel]['sample_rate']),
                                                           DAQmx_Val_Rising, DAQmx_Val_FiniteSamps,
                                                           uInt64(task['sample_num'])))
        else:
            self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle'], "", float64(
                                                            self.settings['analog_input'][channel]['sample_rate']),
                                                           DAQmx_Val_Rising, DAQmx_Val_ContSamps,
                                                           uInt64(task['sample_num'])))

        return task_name

    def setup_DO(self, channels):
        """
        Initializes a arbitrary number of digital output channels to output an arbitrary waveform
        Args:
            channels: List of channels to output, check in self.settings['digital_output'] for available channels
            waveform: 2d array of boolean values to output, with each column giving the output values at a given time
                (the timing given by the sample rate of the channel) with the channels going from top to bottom in
                the column in the order given in channels
            clk_source: the PFI channel of the hardware clock to lock the output to, or "" to use the default
                internal clock

        sets up creates self.DO_taskHandle
        """

        task = {
            'task_handle': None,
            'sample_num': None,
            'sample_rate': None,
            'num_samples_per_channel': None,
            'timeout': None
        }

        task_name = self._add_to_tasklist('do', task)

        if 'digital_output' not in list(self.settings.keys()):
            raise ValueError('This DAQ does not support digital output')
        for c in channels:
            if not c in list(self.settings['digital_output'].keys()):
                raise KeyError('This is not a valid digital output channel')
        task['sample_rate'] = float(
            self.settings['digital_output'][channels[0]]['sample_rate'])  # float prevents truncation in division
        for c in channels:
            if not self.settings['digital_output'][c]['sample_rate'] == task['sample_rate']:
                raise ValueError('All sample rates must be the same')

        lines_list = ''
        for c in channels:
            lines_list += self.settings['device'] + '/port0/line' + str(
                self.settings['digital_output'][c]['channel']) + ','
        lines_list = lines_list[:-1]  # remove the last comma

        self.running = True

        task['task_handle'] = TaskHandle(0)
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle'])))
        self._check_error(self.nidaq.DAQmxCreateDOChan(task['task_handle'],
                                                       lines_list,
                                                       "",
                                                       DAQmx_Val_ChanPerLine))

        return task_name

    def DO_write(self, task_name, output_values):
        task = self.tasklist[task_name]
        sample_num = (numpy.array(output_values).shape)[-1]
        self._check_error(self.nidaq.DAQmxWriteDigitalLines(task['task_handle'],
                                                            int32(sample_num),
                                                            bool32(False),
                                                            float64(-1),
                                                            DAQmx_Val_GroupByChannel,
                                                            np.array(output_values).ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
                                                            None,
                                                            None))

    def read_AI(self, task_name):
        """
        Reads the AI voltage values from the buffer
        Returns: array of ctypes.c_long with the voltage data
        """
        task = self.tasklist[task_name]
        data = (float64 * task['sample_num'])()
        samples_per_channel_read = int32()
        self._check_error(self.nidaq.DAQmxReadAnalogF64(task['task_handle'], task['sample_num'], float64(10.0),
                                                        DAQmx_Val_GroupByChannel, data.ctypes.data,
                                                        task['sample_num'], ctypes.byref(samples_per_channel_read), None))

        return data, samples_per_channel_read


    # run the task specified by task_name
    # todo: AK - should this be threaded? original todo: is this actually blocking? Is the threading actually doing anything? see nidaq cookbook
    def run(self, task_name):
        """
        Runs the task or list of tasks specified in taskname. What 'running' does depends on the type of task that was
        set up, but generally either begins output from a buffer or input to a buffer.

        Args:
            task_name: string identifying task

        """
        #run list of tasks
        if type(task_name) == list:
            for name in task_name:
                task = self.tasklist[name]
                self._check_error(self.nidaq.DAQmxStartTask(task['task_handle']))
        #run single task
        else:
            task = self.tasklist[task_name]
            self._check_error(self.nidaq.DAQmxStartTask(task['task_handle']))

    def waitToFinish(self, task_name):
        """
        Blocks until the task specified by task_name is completed

        Args:
            task_name: string identifying task

        """
        task = self.tasklist[task_name]
        self._check_error(self.nidaq.DAQmxWaitUntilTaskDone(task['task_handle'],
                                                            float64(task['sample_num'] / task['sample_rate'] * 4 + 1)))

    def read(self, task_name):
        if 'ctr' in task_name:
            return(self.read_counter(task_name))
        elif 'ai' in task_name:
            return(self.read_AI(task_name))
        else:
            raise ValueError('This task does not allow reads.')

    def write(self, task_name, output_values):
        if 'do' in task_name:
            self.DO_write(task_name, output_values)
        else:
            raise ValueError('This task does not allow writes.')

    def stop(self, task_name):
        #remove task to be cleared from tasklist
        task = self.tasklist.pop(task_name)

        #special case counters, which create two tasks that need to be cleared
        if 'task_handle_ctr' in list(task.keys()):
            self.nidaq.DAQmxStopTask(task['task_handle_ctr'])
            self.nidaq.DAQmxClearTask(task['task_handle_ctr'])

        self.nidaq.DAQmxStopTask(task['task_handle'])
        self.nidaq.DAQmxClearTask(task['task_handle'])

    def get_analog_voltages(self, channel_list):
        """
        Args:
            channel_list: list (length N) of channels from which to read the voltage, channels are given as strings, e.g. ['ao1', 'ai3']

        Returns:
            list of voltages (length N)

        """
        daq_channels_str = ''
        for channel in channel_list:
            if channel in self.settings['analog_output']:
                daq_channels_str += self.settings['device'] + '/_' + channel + '_vs_aognd, '
            elif (channel in self.settings['analog_input']):
                daq_channels_str += self.settings['device'] + '/' + channel + ', '
        daq_channels_str = daq_channels_str[:-2].encode('ascii')  # strip final comma period
        data = (float64 * len(channel_list))()
        sample_num = 1
        get_voltage_taskHandle = TaskHandle(0)
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(get_voltage_taskHandle)))
        self._check_error(self.nidaq.DAQmxCreateAIVoltageChan(get_voltage_taskHandle, daq_channels_str, "",
                                                              DAQmx_Val_Cfg_Default,
                                                              float64(-10.0), float64(10.0),
                                                              DAQmx_Val_Volts, None))
        self._check_error(self.nidaq.DAQmxReadAnalogF64(get_voltage_taskHandle, int32(sample_num), float64(10.0),
                                                        DAQmx_Val_GroupByChannel, ctypes.byref(data),
                                                        int32(sample_num * len(channel_list)), None, None))
        self._check_error(self.nidaq.DAQmxClearTask(get_voltage_taskHandle))

        for i, channel in enumerate(channel_list):
            # if channel in self.settings['analog_output']:
            data[i] += self.settings['ao_read_offset']

        return [1. * d for d in data]  # return and convert from ctype to python float

    def set_analog_voltages(self, output_dict):
        """

        Args:
            output_dict: dictionary with names of channels as key and voltage as value, e.g. {'ao0': 0.1} or {'0':0.1} for setting channel 0 to 0.1

        Returns: nothing

        """
        # daq API only accepts either one point and one channel or multiple points and multiple channels

        #
        # # make sure the key has the right format, e.g. ao0
        # channels = ['ao'+k.replace('ao','') for k in output_dict.keys()]

        channels = []
        voltages = []
        for k, v in output_dict.items():
            channels.append('ao' + k.replace('ao', ''))  # make sure the key has the right format, e.g. ao0
            voltages.append(v)

        voltages = np.array([voltages]).T
        voltages = (np.repeat(voltages, 2, axis=1))
        # pt = np.transpose(np.column_stack((pt[0],pt[1])))
        # pt = (np.repeat(pt, 2, axis=1))

        task_name = self.setup_AO(channels, voltages)
        self.run(task_name)
        self.waitToFinish(task_name)
        self.stop(task_name)

    def set_digital_output(self, output_dict):
        """

        Args:
            output_dict: dictionary with names of channels as key and voltage as value, e.g. {'do0': True} or {'0':True} for setting channel 0 to True

        Returns: nothing

        """

        channels = []
        values = []
        for k, v in output_dict.items():
            channels.append('do' + k.replace('do', ''))  # make sure the key has the right format, e.g. ao0
            values.append(v)

        print(('channels', channels))
        print(('voltages', values))

        task_name = self.setup_DO(channels)

        self.run(task_name)

        self.DO_write(task_name, values)

        self.stop(task_name)

    def _check_error(self, err):
        """
        Error Checking Routine for DAQmx functions. Pass in the returned values form DAQmx functions (the errors) to get
        an error description. Raises a runtime error
        Args:
            err: 32-it integer error from an NI-DAQmx function

        Returns: a verbose description of the error taken from the nidaq dll

        """
        if err < 0:
            buffer_size = 1000
            buffer = ctypes.create_string_buffer(('\000' * buffer_size).encode('ascii'))
            self.nidaq.DAQmxGetExtendedErrorInfo(ctypes.byref(buffer), buffer_size)
            # raise RuntimeError('nidaq call failed with error %d: %s' % (err, repr(buffer.value)))
            raise RuntimeError('nidaq call failed with error %d: %s' % (err, buffer.value))
        if err > 0:
            buffer_size = 1000
            buffer = ctypes.create_string_buffer(('\000' * buffer_size).encode('ascii'))
            self.nidaq.DAQmxGetErrorString(err, ctypes.byref(buffer), buffer_size)
            raise RuntimeError('nidaq generated warning %d: %s' % (err, repr(buffer.value)))


class NI6259(DAQ):
    """
    This class implements the NI6259 DAQ, which includes 32 AI, 4 AO, and 24 DI/DO channels and inherits basic
    input/output functionality from DAQ. A subset of these channels are accessible here, but more can be added up to
    these limits.
    """
    _DEFAULT_SETTINGS = Parameter([
        Parameter('device', 'Dev1', ['Dev1'], 'Name of DAQ device'),
        Parameter('override_buffer_size', -1, int, 'Buffer size for manual override (unused if -1)'),
        Parameter('ao_read_offset', .005, float, 'Empirically determined offset for reading ao voltages internally'),
        Parameter('analog_output',
                  [
                      Parameter('ao0',
                                [
                                    Parameter('channel', 0, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao1',
                                [
                                    Parameter('channel', 1, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao2',
                                [
                                    Parameter('channel', 2, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao3',
                                [
                                    Parameter('channel', 3, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                )
                  ]
                  ),
        Parameter('analog_input',
                  [
                      Parameter('ai0',
                                [
                                    Parameter('channel', 0, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai1',
                                [
                                    Parameter('channel', 1, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai2',
                                [
                                    Parameter('channel', 2, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai3',
                                [
                                    Parameter('channel', 3, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage')
                                ]
                                ),
                      Parameter('ai4',
                                [
                                    Parameter('channel', 4, list(range(0, 32)), 'input channel'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate'),
                                    Parameter('min_voltage', -10.0, float, 'minimum input voltage'),
                                    Parameter('max_voltage', 10.0, float, 'maximum input voltage (V)')
                                ]
                                )
                  ]
                  ),
        Parameter('digital_input',
                  [
                      Parameter('ctr0',
                                [
                                    Parameter('input_channel', 0, list(range(0, 32)), 'channel for counter signal input'),
                                    Parameter('counter_PFI_channel', 8, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('gate_PFI_channel', 14, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('clock_PFI_channel', 13, list(range(0, 32)), 'PFI for clock channel output'),
                                    Parameter('clock_counter_channel', 1, [0, 1], 'channel for clock output'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)')
                                ]
                                ),
                      Parameter('ctr1',
                                [
                                    Parameter('input_channel', 1, list(range(0, 32)), 'channel for counter signal input'),
                                    Parameter('counter_PFI_channel', 3, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('gate_PFI_channel', 14, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('clock_PFI_channel', 12, list(range(0, 32)), 'PFI for clock channel output'),
                                    Parameter('clock_counter_channel', 0, [0, 1], 'channel for clock output'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)')
                                ]
                                )
                  ]
                  ),
        Parameter('digital_output',
                  [
                      Parameter('do0',
                                [
                                    Parameter('channel', 0, list(range(0, 16)), 'channel'),
                                    # Parameter('value', False, bool, 'value')
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)')
                                    # Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    # Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('do8',
                                [
                                    Parameter('channel', 8, list(range(0, 16)), 'channel'),
                                    # Parameter('value', False, bool, 'value')
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)')
                                    # Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    # Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                )
                  ]
                  )
    ])

class NI9263(DAQ):
    """
    This class implements the NI9263 DAQ, which includes 4 AO channels. It inherits output functionality from the DAQ
    class.
    """
    _DEFAULT_SETTINGS = Parameter([
        Parameter('device', 'cDAQ1Mod1', ['cDAQ9184-1BA7633Mod3', 'cDAQ9184-1BA7633Mod4', 'cDAQ9184-1BA7633Mod1', 'cDAQ1Mod1'],
                  'Name of DAQ device - check in NiMax'),
        Parameter('override_buffer_size', -1, int, 'Buffer size for manual override (unused if -1)'),
        Parameter('ao_read_offset', .005, float, 'Empirically determined offset for reading ao voltages internally'),
        Parameter('analog_output',
                  [
                      Parameter('ao0',
                                [
                                    Parameter('channel', 0, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao1',
                                [
                                    Parameter('channel', 1, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao2',
                                [
                                    Parameter('channel', 2, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                ),
                      Parameter('ao3',
                                [
                                    Parameter('channel', 3, [0, 1, 2, 3], 'output channel'),
                                    Parameter('sample_rate', 1000.0, float, 'output sample rate (Hz)'),
                                    Parameter('min_voltage', -10.0, float, 'minimum output voltage (V)'),
                                    Parameter('max_voltage', 10.0, float, 'maximum output voltage (V)')
                                ]
                                )
                  ]
                  )
    ])

class NI9402(DAQ):
    """
    This class implements the NI9263 DAQ, which includes 4 AO channels. It inherits output functionality from the DAQ
    class.
    """
    _DEFAULT_SETTINGS = Parameter([
        Parameter('device', 'cDAQ1', ['cDAQ1', 'cDAQ9184-1BA7633', 'cDAQ9188-1BFB6F2'], 'Name of DAQ device - check in NiMax'),
        Parameter('module', 'Mod2', ['Mod1', 'Mod2', 'Mod3', 'Mod4', 'Mod5', 'Mod6', 'Mod7', 'Mod8']),
        Parameter('override_buffer_size', -1, int, 'Buffer size for manual override (unused if -1)'),
        Parameter('ao_read_offset', .005, float, 'Empirically determined offset for reading ao voltages internally'),
        Parameter('digital_input',
                  [
                      Parameter('ctr0',
                                [
                                    Parameter('input_channel', 0, list(range(0, 32)), 'channel for counter signal input'),
                                    Parameter('counter_PFI_channel', 0, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('gate_PFI_channel', 3, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('clock_PFI_channel', 1, list(range(0, 32)), 'PFI for clock channel output'),
                                    Parameter('clock_counter_channel', 2, list(range(0, 32)), 'channel for clock output'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)')
                                ]
                                ),
                      Parameter('ctr2',
                                [
                                    Parameter('input_channel', 2, list(range(0, 32)), 'channel for counter signal input'),
                                    Parameter('counter_PFI_channel', 1, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('gate_PFI_channel', 0, list(range(0, 32)), 'PFI for counter channel input'),
                                    Parameter('clock_PFI_channel', 2, list(range(0, 32)), 'PFI for clock channel output'),
                                    Parameter('clock_counter_channel', 3, list(range(0, 32)), 'channel for clock output'),
                                    Parameter('sample_rate', 1000.0, float, 'input sample rate (Hz)')
                                ]
                                )
                  ]
                  ),
    ])

    def setup_counter(self, channel, sample_num, continuous_acquisition=False):
        """
        Initializes a hardware-timed digital counter, bound to a hardware clock
        Args:
            channel: digital channel to initialize for read in
            sample_num: number of samples to read in for finite operation, or number of samples between
                       reads for continuous operation (to set buffer size)
            continuous_acquisition: run in continuous acquisition mode (ex for a continuous counter) or
                                    finite acquisition mode (ex for a scan, where the number of samples needed
                                    is known a priori)

        Returns: source of clock that this method sets up, which can be given to another function to synch that
        input or output to the same clock

        """

        # Note that for this counter, we have two tasks. The normal 'task_handle' corresponds to the clock, and this
        # is the task which is started when run is called. The second 'task_handle_ctr' corresponds to the counter,
        # and this waits for the clock and will be started simultaneously.
        task = {
            'task_handle': None,
            'task_handle_ctr': None,
            'counter_out_PFI_str': None,
            'sample_num': None,
            'sample_rate': None,
            'num_samples_per_channel': None,
            'timeout': None
        }

        task_name = self._add_to_tasklist('ctr', task)

        if 'digital_input' not in list(self.settings.keys()):
            raise ValueError('This DAQ does not support digital input')
        if not channel in list(self.settings['digital_input'].keys()):
            raise KeyError('This is not a valid digital input channel')

        channel_settings = self.settings['digital_input'][channel]
        self.running = True
        task['sample_num'] = sample_num
        task['sample_rate'] = float(channel_settings['sample_rate'])
        if not continuous_acquisition:
            task['num_samples_per_channel'] = task['sample_num']
        else:
            task['num_samples_per_channel'] = -1
        task['timeout'] = float64(5 * (1 / task['sample_rate']) * task['sample_num'])
        input_channel_str = (self.settings['device'] + self.settings['module'] + '/' + channel).encode('ascii')
        task['counter_out_PFI_str'] = ('/' + self.settings['device'] + '/Ctr' + str(channel_settings['clock_counter_channel']) + 'InternalOutput').encode('ascii')
        counter_out_str = (self.settings['device'] + self.settings['module'] + '/ctr' + str(channel_settings['clock_counter_channel'])).encode('ascii')
        task['task_handle_ctr'] = TaskHandle(0)
        task['task_handle'] = TaskHandle(1)

        # set up clock
        self._dig_pulse_train_cont(task, .5, counter_out_str)
        # set up counter using clock as reference
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle_ctr'])))
        self._check_error(self.nidaq.DAQmxCreateCICountEdgesChan(task['task_handle_ctr'],
                                                                 input_channel_str, "", DAQmx_Val_Rising, 0,
                                                                 DAQmx_Val_CountUp))

        if not continuous_acquisition:
            self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle_ctr'], task['counter_out_PFI_str'],
                                                               float64(task['sample_rate']), DAQmx_Val_Rising,
                                                               DAQmx_Val_FiniteSamps, uInt64(task['sample_num'])))
        else:
            self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle_ctr'], task['counter_out_PFI_str'],
                                                               float64(task['sample_rate']), DAQmx_Val_Rising,
                                                               DAQmx_Val_ContSamps, uInt64(task['sample_num'])))
        # if (self.settings['override_buffer_size'] > 0):
        #     self._check_error(self.nidaq.DAQmxCfgInputBuffer(self.DI_taskHandleCtr, uInt64(self.settings['override_buffer_size'])))
        # self._check_error(self.nidaq.DAQmxCfgInputBuffer(self.DI_taskHandleCtr, uInt64(sampleNum)))

        self._check_error(self.nidaq.DAQmxStartTask(task['task_handle_ctr']))


        return task_name

if __name__ == '__main__':
    pass
    # daq, failed = Instrument.load_and_append({'daq': NI9263, 'daq_in': NI6259})
    # print('FAILED', failed)
    # print(daq['daq'].settings)
    #
    # daq['daq'].device = 'cDAQ9184-1BA7633Mod3'
    # print('------ daq ------')
    # print(daq['daq'])
    # print('------ daq_in ------')
    # print(daq['daq_in'])
    #
    # print('------')
    #
    # vout = -1
    # daq['daq'].set_analog_voltages({'ao0': vout})
    # print('SET', vout)
    #
    # print(daq['daq_in'])
    # print('GET', daq['daq_in'].get_analog_voltages(['ai1', 'ai2']))

    # daq, failed = Instrument.load_and_append({'daq_in': NI6259})
    # daq['daq_in'].set_digital_output({'do0': False})
