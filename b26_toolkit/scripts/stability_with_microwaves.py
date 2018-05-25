from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.instruments import MicrowaveGenerator, NI6259
from collections import deque
import time
# from b26_toolkit.pylabcontrol.plotting.plots_1d import plot_esr
# from b26_toolkit.pylabcontrol.data_processing.esr_signal_processing import fit_esr

class Stability_With_Microwaves(Script):
    """
    This script alternatingly outputs two microwave frequencies and measured fluourescent counts for each. The difference in counts is output continuously.
    This serves as a real time signal for the magnetic field splitting if one frequency is near the dip of the esr spectrum and the other far away as a reference.
    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -45.0, float, 'output power (dBm)'),
        Parameter('freq', 2.82e9, float, 'first frequency (Hz)'),
        Parameter('measurement_time', 0.5, float, 'measurement time of fluorescent counts for each of the three windows (s)'),
        Parameter('settle_time', 10, float, 'dead time after switching frequency (ms)'),
        Parameter('interval_time', .01, float, 'binning time for each measurement window (s)')
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

        daq = self.instruments['daq']['instance']
        mw_gen = self.instruments['microwave_generator']['instance']

        freq = self.settings['freq']
        measurement_time = self.settings['measurement_time']
        interval_time = self.settings['interval_time']
        settle_time = self.settings['settle_time']

        sample_rate = float(1) / interval_time
        self.instruments['daq']['instance'].settings['digital_input']['ctr0']['sample_rate'] = sample_rate

        self.instruments['microwave_generator']['instance'].update({'amplitude': self.settings['power_out']})
        self.instruments['microwave_generator']['instance'].update({'enable_modulation': False})

        self.data = {'time': [], 'counts_before': [], 'counts_during': [], 'counts_after': []}

        for window in ['counts_before', 'counts_during', 'counts_after']:

            if window == 'counts_during':
                mw_gen.update({'enable_output': True})
            else:
                mw_gen.update({'enable_output': False})

            # create DAQ task
            task = daq.setup_counter("ctr0", sample_num=int(measurement_time/interval_time) + 1)  # adds +1 so diff returns proper length array
            # run daq task to read APD
            daq.run(task)
            raw_data, _ = daq.read_counter(task)  # blocking read
            daq.stop(task)

            self.data['time'] = np.linspace(0,measurement_time, int(measurement_time/interval_time))
            self.data[window] = np.diff(raw_data) * sample_rate / 1000


        self.progress = 50
        self.updateProgress.emit(self.progress)


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
        print(('time', self.data['time']))
        print(('counts_before', self.data['counts_before']))

        axes_list[0].plot(self.data['time'], self.data['counts_before'], label = 'before_mw')
        axes_list[0].plot(self.data['time'], self.data['counts_during'], label = 'during_mw')
        axes_list[0].plot(self.data['time'], self.data['counts_after'], label = 'after_mw')
        axes_list[0].legend()
        axes_list[0].set_xlabel('time (s)')
        axes_list[0].set_ylabel('signal (kCounts/s)')