import ctypes
import numpy as np
import datetime

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

class DAQ():
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

    dll_path = 'C:/Windows/System32/nicaiu.dll'
    nidaq = ctypes.WinDLL(dll_path)  # load the DLL

    tasklist = {}
    tasknum = 0

    settings = {}

    def _add_to_tasklist(self, name, task):
        """
        Adds task to a list of currently initialized tasks
        Args:
            name: name defining the type of task
            task: dictionary containing the task information

        Returns: name of task added to tasklist

        """
        matching = [x for x in self.tasklist if name in x]
        if not matching:
            task_name = name + '000'
        else:
            last_task = sorted(matching)[-1]
            task_name = name + '{0:03d}'.format(int(last_task[-3:])+1)
        self.tasklist.update({task_name: task})
        return task_name

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
        task['num_samples_per_channel'] = task['sample_num']

        task['timeout'] = float64(5 * (1 / task['sample_rate']) * task['sample_num'])
        input_channel_str = self.settings['device'] + self.settings['module'] + '/' + channel
        task['counter_out_PFI_str'] = '/' + self.settings['device'] + '/Ctr' + str(channel_settings['clock_counter_channel']) + 'InternalOutput'
        counter_out_str = self.settings['device'] + self.settings['module'] + '/' + 'ctr' + str(channel_settings['clock_counter_channel'])
        task['task_handle_ctr'] = TaskHandle(0)
        task['task_handle'] = TaskHandle(1)

        # set up clock
        self._dig_pulse_train_cont(task, .5, counter_out_str)
        # set up counter using clock as reference
        self._check_error(self.nidaq.DAQmxCreateTask("", ctypes.byref(task['task_handle_ctr'])))
        self._check_error(self.nidaq.DAQmxCreateCICountEdgesChan(task['task_handle_ctr'],
                                                                 input_channel_str, "", DAQmx_Val_Rising, 0,
                                                                 DAQmx_Val_CountUp))

        self._check_error(self.nidaq.DAQmxCfgSampClkTiming(task['task_handle_ctr'], task['counter_out_PFI_str'],
                                                           float64(task['sample_rate']), DAQmx_Val_Rising,
                                                           DAQmx_Val_FiniteSamps, uInt64(task['sample_num'])))
        # self._check_error(self.nidaq.DAQmxTaskControl(task['task_handle_ctr'], 3))


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
        # self._check_error(self.nidaq.DAQmxTaskControl(task['task_handle'], 3))

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
        channel_list = ''
        for c in channels:
            channel_list += self.settings['device'] + '/' + c + ','
        channel_list = channel_list[:-1]
        self.running = True
        # special case 1D waveform since length(waveform[0]) is undefined
        if (len(np.shape(waveform)) == 2):
            numChannels = len(waveform)
            task['sample_num'] = len(waveform[0])
        else:
            task['sample_num'] = len(waveform)
            numChannels = 1
        task['task_handle'] = TaskHandle(self.tasknum)
        self.tasknum += 1
        # special case 1D waveform since length(waveform[0]) is undefined
        # converts python array to ctypes array
        if (len(np.shape(waveform)) == 2):
            data = np.zeros((numChannels, task['sample_num']),dtype=np.float64)
            for i in range(numChannels):
                for j in range(task['sample_num']):
                    data[i, j] = waveform[i, j]
        else:
            data = np.zeros((task['sample_num']), dtype=np.float64)
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
        # self._check_error(self.nidaq.DAQmxTaskControl(task['task_handle'], 3))

        return task_name

    def run(self, task_name):
        """
        Runs the task or list of tasks specified in taskname. What 'running' does depends on the type of task that was
        set up, but generally either begins output from a buffer or input to a buffer.

        Args:
            task_name: string identifying task

        """
        #run list of tasks
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
        return(self.read_counter(task_name))

    def stop(self, task_name):
        #remove task to be cleared from tasklist
        # task = self.tasklist.pop(task_name)
        task = self.tasklist[task_name]

        #special case counters, which create two tasks that need to be cleared
        if 'task_handle_ctr' in list(task.keys()):
            self.nidaq.DAQmxStopTask(task['task_handle_ctr'])
            self.nidaq.DAQmxClearTask(task['task_handle_ctr'])

        self.nidaq.DAQmxStopTask(task['task_handle'])
        self.nidaq.DAQmxClearTask(task['task_handle'])

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
            buffer = ctypes.create_string_buffer('\000' * buffer_size)
            # self.nidaq.DAQmxGetErrorString(err, ctypes.byref(buffer), buffer_size)
            self.nidaq.DAQmxGetExtendedErrorInfo(ctypes.byref(buffer), buffer_size)
            # raise RuntimeError('nidaq call failed with error %d: %s' % (err, repr(buffer.value)))
            raise RuntimeError('nidaq call failed with error %d: %s' % (err, buffer.value))
        if err > 0:
            buffer_size = 1000
            buffer = ctypes.create_string_buffer('\000' * buffer_size)
            self.nidaq.DAQmxGetErrorString(err, ctypes.byref(buffer), buffer_size)
            raise RuntimeError('nidaq generated warning %d: %s' % (err, repr(buffer.value)))

class NI9263(DAQ):
    """
    This class implements the NI9263 DAQ, which includes 4 AO channels. It inherits output functionality from the DAQ
    class.
    """
    def __init__(self):
        self.settings = {'device': 'cDAQ9184-1BA7633Mod1', #'cDAQ9188-1BFB6F2Mod1',
                    'override_buffer_size': -1,
                    'ao_read_offset': .005,
                    'analog_output': {'ao0': {'channel': 0, 'sample_rate': 1000.0, 'min_voltage': -10.0, 'max_voltage': 10.0},
                                    'ao1': {'channel': 1, 'sample_rate': 1000.0, 'min_voltage': -10.0, 'max_voltage': 10.0}}
        }

class NI9402(DAQ):
    """
    This class implements the NI9402 DAQ, which uses 4 counter channels from the cDAQ chassis. It inherits output
    functionality from the DAQ class.
    """
    def __init__(self):
        self.settings = {'device': 'cDAQ9188-1BFB6F2',
                    'module': 'Mod2',
                    'override_buffer_size': -1,
                    'ao_read_offset': .005,
                    'digital_input': {'ctr0': {'input_channel': 0, 'counter_PFI_channel': 0, 'clock_PFI_channel': 1, 'clock_counter_channel': 2, 'sample_rate': 1000.0}}
                    }

_INSTRUMENTS = {'daq_out': NI9263, 'daq_in': NI9402}

def pts_to_extent(pta, ptb, roi_mode):
    """

    Args:
        pta: point a
        ptb: point b
        roi_mode:   mode how to calculate region of interest
                    corner: pta and ptb are diagonal corners of rectangle.
                    center: pta is center and ptb is extend or rectangle

    Returns: extend of region of interest [xVmin, xVmax, yVmax, yVmin]

    """
    if roi_mode == 'corner':
        xVmin = min(pta['x'], ptb['x'])
        xVmax = max(pta['x'], ptb['x'])
        yVmin = min(pta['y'], ptb['y'])
        yVmax = max(pta['y'], ptb['y'])
    elif roi_mode == 'center':
        xVmin = pta['x'] - float(ptb['x']) / 2.
        xVmax = pta['x'] + float(ptb['x']) / 2.
        yVmin = pta['y'] - float(ptb['y']) / 2.
        yVmax = pta['y'] + float(ptb['y']) / 2.
    return [xVmin, xVmax, yVmax, yVmin]

def RUN_ME():
    # start time
    print((datetime.datetime.now()))

    instruments = {'daq_out': NI9263(), 'daq_in': NI9402()}

    settings = {'point_a': {'x': 0, 'y': 0},
                'point_b': {'x': 1, 'y': 1},
                'RoI_mode': 'center',
                'num_points': {'x': 126, 'y': 126},
                'time_per_pt': .002,
                'settle_time': .0002,
                'max_counts_plot': -1,
                'DAQ_channels': {'x_ao_channel': 'ao0', 'y_ao_channel': 'ao1', 'counter_channel': 'ctr0'}
                }

    #initialize scan
    clockAdjust = int(
        (settings['time_per_pt'] + settings['settle_time']) / settings['settle_time'])

    [xVmin, xVmax, yVmax, yVmin] = pts_to_extent(settings['point_a'], settings['point_b'], settings['RoI_mode'])

    x_array = np.repeat(
        np.linspace(xVmin, xVmax, settings['num_points']['x'], endpoint=True),
        clockAdjust)
    y_array = np.linspace(yVmin, yVmax, settings['num_points']['y'], endpoint=True)
    sample_rate = float(1) / settings['settle_time']
    instruments['daq_out'].settings['analog_output'][
        settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
    instruments['daq_out'].settings['analog_output'][
        settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate
    instruments['daq_in'].settings['digital_input'][
        settings['DAQ_channels']['counter_channel']]['sample_rate'] = sample_rate
    data = {'image_data': np.zeros((settings['num_points']['y'], settings['num_points']['x'])),
                 'bounds': [xVmin, xVmax, yVmin, yVmax]}

    data['extent'] = [xVmin, xVmax, yVmax, yVmin]


    #begin scan
    for yNum in range(0, len(y_array)):

        # set galvo to initial point of next line
        initPt = [x_array[0], y_array[yNum]]
        instruments['daq_out'].set_analog_voltages(
            {settings['DAQ_channels']['x_ao_channel']: initPt[0],
             settings['DAQ_channels']['y_ao_channel']: initPt[1]})

        # initialize APD thread
        ctrtask = instruments['daq_in'].setup_counter(
            settings['DAQ_channels']['counter_channel'],
            len(x_array) + 1)
        aotask = instruments['daq_out'].setup_AO([settings['DAQ_channels']['x_ao_channel']],
                                                              x_array, ctrtask)


        # start counter and scanning sequence
        instruments['daq_out'].run(aotask)
        instruments['daq_in'].run(ctrtask)
        instruments['daq_out'].waitToFinish(aotask)
        instruments['daq_out'].stop(aotask)
        xLineData, _ = instruments['daq_in'].read(ctrtask)
        instruments['daq_in'].stop(ctrtask)
        diffData = np.diff(xLineData)

        summedData = np.zeros(len(x_array) / clockAdjust)
        for i in range(0, int((len(x_array) / clockAdjust))):
            summedData[i] = np.sum(
                diffData[(i * clockAdjust + 1):(i * clockAdjust + clockAdjust - 1)])
        # also normalizing to kcounts/sec
        data['image_data'][yNum] = summedData * (.001 / settings['time_per_pt'])

        return

    # end time
    print((datetime.datetime.now()))

def RUN_ME_output_only():
    # start time
    print((datetime.datetime.now()))

    instruments = {'daq_out': NI9263()}

    settings = {'point_a': {'x': 0, 'y': 0},
                'point_b': {'x': 1, 'y': 1},
                'RoI_mode': 'center',
                'num_points': {'x': 126, 'y': 126},
                'time_per_pt': .002,
                'settle_time': .0002,
                'max_counts_plot': -1,
                'DAQ_channels': {'x_ao_channel': 'ao0', 'y_ao_channel': 'ao1', 'counter_channel': 'ctr0'}
                }

    #initialize scan
    clockAdjust = int(
        (settings['time_per_pt'] + settings['settle_time']) / settings['settle_time'])

    [xVmin, xVmax, yVmax, yVmin] = pts_to_extent(settings['point_a'], settings['point_b'], settings['RoI_mode'])

    x_array = np.repeat(
        np.linspace(xVmin, xVmax, settings['num_points']['x'], endpoint=True),
        clockAdjust)
    y_array = np.linspace(yVmin, yVmax, settings['num_points']['y'], endpoint=True)
    sample_rate = float(1) / settings['settle_time']
    instruments['daq_out'].settings['analog_output'][
        settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = sample_rate
    instruments['daq_out'].settings['analog_output'][
        settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = sample_rate

    # initialize APD thread
    # aotask = instruments['daq_out'].setup_AO([settings['DAQ_channels']['x_ao_channel']],
    #                                          x_array)

    #begin scan
    for yNum in range(0, len(y_array)):

        # set galvo to initial point of next line
        initPt = [x_array[0], y_array[yNum]]
        instruments['daq_out'].set_analog_voltages(
            {settings['DAQ_channels']['x_ao_channel']: initPt[0],
             settings['DAQ_channels']['y_ao_channel']: initPt[1]})

        aotask = instruments['daq_out'].setup_AO([settings['DAQ_channels']['x_ao_channel']],
                                                 x_array)

        # start counter and scanning sequence
        instruments['daq_out'].run(aotask)
        instruments['daq_out'].waitToFinish(aotask)
        instruments['daq_out'].stop(aotask)

        # return

    # end time
    print((datetime.datetime.now()))

RUN_ME_output_only()