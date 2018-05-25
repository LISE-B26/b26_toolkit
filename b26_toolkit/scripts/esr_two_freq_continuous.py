from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.instruments import MicrowaveGenerator, NI6259
from collections import deque
import time
# from b26_toolkit.pylabcontrol.plotting.plots_1d import plot_esr
# from b26_toolkit.pylabcontrol.data_processing.esr_signal_processing import fit_esr

class ESRTwoFreqContinuous(Script):
    """
    This script alternatingly outputs two microwave frequencies and measured fluourescent counts for each. The difference in counts is output continuously.
    This serves as a real time signal for the magnetic field splitting if one frequency is near the dip of the esr spectrum and the other far away as a reference.
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('freq_1', 2.82e9, float, 'first frequency (Hz)'),
        Parameter('freq_2', 2.92e9, float, 'second frequency (Hz)'),
        Parameter('measurement_time', 0.01, float, 'measurement time of fluorescent counts (s)'),
        Parameter('settle_time', 10, float, 'dead time after switching frequency (ms)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('max_points', 100, int, 'number of points to display if 0 show all'),
        Parameter('range_type', 'freq_1_2', ['freq_1_2', 'freq_1_delta'],
                  'freq_1_2: measure at frequency 1 and frequency 2 freq_0_delta: measure at frequency 1 and (frequency 1 - frequency 2)'),
    ]

    _INSTRUMENTS = {
        'microwave_generator': MicrowaveGenerator,
        'daq': NI6259
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts = None, name=None, settings=None, log_function=None, data_path = None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path = data_path)
    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        def set_freq_and_read_daq(freq, settle_time, daq, mw_gen):
            """
            reads a frequency section from the DAQ

            Args:
                freq: output frequency of microwave_generator
                settle_time:  time it takes to set the frequency in ms, typically 10ms
                daq: instance of the daq, e.g. self.instruments['daq']['instance']
                mw_gen: instance of the mw generator, e.g. self.instruments['microwave_generator']['instance']

            Returns: count rate from daq

            """

            # create DAQ task
            task = daq.setup_counter("ctr0", sample_num=10) #comment here
            # set frequency
            mw_gen.update({'frequency': float(freq)})
            # wait for instrument to apply settings
            time.sleep(1e-3*settle_time)
            # run daq task to read APD
            daq.run(task)
            raw_data, _ = daq.read_counter(task) #blocking read
            daq.stop(task)

            count_rate = np.mean(np.diff(raw_data)) * sample_rate / 1000
            return count_rate


        measurement_time = self.settings['measurement_time']
        settle_time = self.settings['settle_time']


        if self.settings['range_type'] == 'freq_1_2':
            freq_1 = self.settings['freq_1']
            freq_2 = self.settings['freq_2']
        elif self.settings['range_type'] == 'freq_1_delta':
            freq_1 = self.settings['freq_1']
            freq_2 = self.settings['freq_1'] - self.settings['freq_2']
        else:
            raise NotImplementedError('unknown setting for range_type')


        sample_rate = float(10) / measurement_time
        self.instruments['daq']['instance'].settings['digital_input']['ctr0']['sample_rate'] = sample_rate

        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        # self.instruments['microwave_generator']['instance'].update({'modulation_type': 'FM'})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})

        self.progress = 50

        if self.settings['max_points']>0:
            esr_data = [deque(maxlen = self.settings['max_points']), deque(maxlen = self.settings['max_points'])]
        else:
            esr_data = [deque(), deque()]
        self.data = {'frequency': [], 'data': esr_data, 'fit_params': []}


        while self._abort is False:
            count_rate1 = set_freq_and_read_daq(freq_1, settle_time, self.instruments['daq']['instance'], self.instruments['microwave_generator']['instance'])
            count_rate2 = set_freq_and_read_daq(freq_2, settle_time, self.instruments['daq']['instance'], self.instruments['microwave_generator']['instance'])

            esr_data[0].append(count_rate1)
            esr_data[1].append(count_rate2)


            self.updateProgress.emit(self.progress)


        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})


    def _calc_progress(self):
        #COMMENT_ME
        self.progress = 50
        return int(self.progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """
        if data is None:
            data = np.array(self.data['data'])

        max_points = self.settings['max_points']

        contrast = 100.*(data[1] - data[0])/(data[0] + data[1])

        axes_list[0].plot(contrast, 'b')
        axes_list[0].plot([0,max_points],[0,0],'k-')
        axes_list[0].set_title('contrast')
        axes_list[1].plot(data[0], 'r')
        axes_list[1].plot(data[1], 'b')
        axes_list[1].set_title('esr from each freq. (1=r), (2=b)')
        axes_list[0].set_xlabel('time (arb units)')
        axes_list[1].set_xlabel('time (arb units)')
        axes_list[1].set_ylabel('kCounts/s')
        axes_list[0].set_ylabel('contrast (%)')


    def _update_plot(self, axes_list, data = None):
        if data is None:
            data = np.array(self.data['data'])

        contrast = 100.*(data[1] - data[0])/(data[0] + data[1])

        axes_list[0].lines[0].set_xdata(list(range(0, len(data[0]))))
        axes_list[0].lines[0].set_ydata(contrast)

        axes_list[0].relim()
        axes_list[0].autoscale_view()

        axes_list[1].lines[0].set_xdata(list(range(0,len(data[0]))))
        axes_list[1].lines[0].set_ydata(data[0])
        axes_list[1].lines[1].set_xdata(list(range(0, len(data[1]))))
        axes_list[1].lines[1].set_ydata(data[1])

        axes_list[0].set_xlabel('time (arb units)')
        axes_list[1].set_xlabel('time (arb units)')
        axes_list[1].set_ylabel('kCounts/s')
        axes_list[0].set_ylabel('contrast (%)')

        axes_list[1].relim()
        axes_list[1].autoscale_view()