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
        task['task_handle'] = TaskHandle(0)
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

    def stop(self, task_name):
        #remove task to be cleared from tasklist
        task = self.tasklist.pop(task_name)

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
        return
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
        self.settings = {'device': 'cDAQ9188-1BFB6F2Mod1',
                    'analog_output': {'ao0': {'channel': 0, 'sample_rate': 1000.0, 'min_voltage': -10.0, 'max_voltage': 10.0},
                                    'ao1': {'channel': 1, 'sample_rate': 1000.0, 'min_voltage': -10.0, 'max_voltage': 10.0}}
        }

def RUN_ME():
    instruments = {'daq_out': NI9263()}

    settings = {'DAQ_channels': {'x_ao_channel': 'ao0', 'y_ao_channel': 'ao1'}}

    #begin scan
    # start time
    print((datetime.datetime.now()))

    # set ao0 and ao1 to 0 V
    initPt = [0, 0]
    instruments['daq_out'].set_analog_voltages(
        {settings['DAQ_channels']['x_ao_channel']: initPt[0],
         settings['DAQ_channels']['y_ao_channel']: initPt[1]})

    # end time
    print((datetime.datetime.now()))

RUN_ME()