from PyLabControl.src.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.src.instruments import MicrowaveGenerator, NI6259
from collections import deque
import time
# from b26_toolkit.src.plotting.plots_1d import plot_esr
# from b26_toolkit.src.data_processing.esr_signal_processing import fit_esr

class ESRTwoFreqContinuous(Script):
    """
    This script alternatingly outputs two microwave frequencies and measured fluourescent counts for each. The difference in counts is output continuously.
    This serves as a real time signal for the magnetic field splitting if one frequency is near the dip of the esr spectrum and the other far away as a reference.
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('esr_avg', 50, int, 'number of esr averages'),
        Parameter('freq_1', 2.82e9, float, 'first frequency (Hz)'),
        Parameter('freq_2', 2.92e9, float, 'second frequency (Hz)'),
        Parameter('measurement_time', 0.01, float, 'measurement time of fluorescent counts (s)'),
        Parameter('settle_time', 10, float, 'dead time after switching frequency (ms)'),
        Parameter('turn_off_after', False, bool, 'if true MW output is turned off after the measurement'),
        Parameter('max_points', 100, int, 'number of points to display if 0 show all')
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
            task = daq.setup_counter("ctr0", sample_num=1)
            # set frequency
            mw_gen.update({'frequency': float(freq)})
            # wait for instrument to apply settings
            time.sleep(1e-3*settle_time)
            # run daq task to read APD
            daq.run(task)
            daq.waitToFinish(task)
            daq.stop(task)

            raw_data, _ = daq.read_counter(task)
            count_rate = raw_data * sample_rate
            return count_rate


        self.lines = []

        take_ref = self.settings['take_ref']
        measurement_time = self.settings['measurement_time']
        settle_time = self.settings['settle_time']

        freq_1 = self.settings['freq_1']
        freq_2 = self.settings['freq_2']


        sample_rate = float(1) / measurement_time
        normalization = self.settings['integration_time']/.001 # to convert from ms to s?
        self.instruments['daq']['instance'].settings['digital_input'][self.settings['counter_channel']]['sample_rate'] = sample_rate

        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'modulation_type': 'FM'})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': True})

        task = self.instruments['daq']['instance'].setup_counter("ctr0", sample_num = 1)


        esr_data = deque()
        if self.settings['max_points']>0:
            esr_data.maxlen = self.settings['max_points']
        self.data = {'frequency': [], 'data': esr_data, 'fit_params': []}


        while self._abort is False:
            count_rate1 = set_freq_and_read_daq(freq_1, settle_time, self.instruments['daq']['instance'], self.instruments['microwave_generator']['instance'])
            count_rate2 = set_freq_and_read_daq(freq_2, settle_time, self.instruments['daq']['instance'], self.instruments['microwave_generator']['instance'])

            esr_data.append([count_rate1, count_rate2])




        if self.settings['turn_off_after']:
            self.instruments['microwave_generator']['instance'].update({'enable_output': False})


    def _calc_progress(self, scan_num):
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

        axes_list[0].plot(data[0:])
        axes_list[1].plot(data[1:])
