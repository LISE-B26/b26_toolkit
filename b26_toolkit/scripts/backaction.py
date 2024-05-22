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

from collections import deque
import numpy as np
import datetime
import time

from b26_toolkit.plotting.plots_1d import plot_psd
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.hf2li_scripts import LockInDaqRead
from b26_toolkit.scripts.pulse_sequences.pi_pulse_train import PiPulseTrainBackaction
from b26_toolkit.instruments import MicrowaveGenerator, Hf2Li
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace
from scipy.signal import periodogram as periodogram

class BackactionGeneric(Script):
    """
    Runs a pulse sequence while varying some parameter, and measures the mechanical response using the lock-in DAQ
    """
    _DEFAULT_SETTINGS = [
        Parameter('param_start', 2.87e9, float, 'start value of parameter sweep'),
        Parameter('param_stop', 1e8, float, 'end value of parameter sweep'),
        Parameter('param_points', 50, int, 'number of parameter points in scan'),
        Parameter('param_switching_time', .01, float, 'time wait after changing parameter value (s)'),
        Parameter('num_averages', 1, int, 'number of sweeps to perform and average over')
    ]

    _INSTRUMENTS = {'mw_gen': MicrowaveGenerator, 'hf2li': Hf2Li}

    _SCRIPTS = {'backaction_pulse_seq': PiPulseTrainBackaction, 'lock_in_daq_read': LockInDaqRead}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments, log_function=log_function, data_path=data_path)
        self.params = []
        self.sweep_params = {}

    def define_sweep_parameters(self):
        """
        Define name of sweep parameters. Since this script is generic, it is coded with variables like 'param_start'.
        This function redirects the code to look for the corresponding settings in _DEFAULT_SETTINGS
        :return:
        """
        self.sweep_params['param_start'] = self.settings['param_start']
        self.sweep_params['param_stop'] = self.settings['param_stop']
        self.sweep_params['param_points'] = self.settings['param_points']
        self.sweep_params['param_switching_time'] = self.settings['param_switching_time']

    def _configure_instruments_start_of_script(self):
        """
        Configure instruments at the beginning of a script, e.g. enable IQ modulation
        :return: None
        """
        pass

    def _configure_instruments_end_of_script(self):
        """
        Configure instruments right before the script finishes, e.g. turn off function generators
        :return: None
        """
        pass

    def _configure_instruments_start_of_sweep(self, param_current):
        """
        Configure instruments before running a sweep. For example, change MW frequency
        :return: None
        """
        raise NotImplementedError

    def _configure_subscripts_start_of_sweep(self, param_current):
        """

        :param param_current:
        :return:
        """

    def _configure_param_array(self):
        """
        Contruct the parameter array and store it in a variable called 'mw_frequencies'. Despite the naming, it's just a list of parameters to be swept;
        it can be a MHz frequency array for a function generator, or even some constant voltage to a LED diode
        :return:
        """

        raise NotImplementedError

    def _function(self, in_data=None):
        """
        The core function can be configured in two ways. In the first, it iterates over each MW freq/param and calls the _run_sweep() function
        from pulsed_experiment_base_script. In this _run_sweep() function, we tell the PB to run the same pulse sequence for one averaging blcok
        and then read out the counts (it's a "sweep" over a list of one pulse seq, as opposed to the usual sweep over many pulse seqs with different
        pulse durations (tau).

        In the second way, we program and start the PB with averaging block size * num of freq points sequences, and read the DAQ counts AFTER iterating
        over all the freq values. This eliminates a lot of the overhead from having to program/read the DAQ repeatedly for each freq point, at the cost
        of less frequent updates. To synchronize the pulse seqs with the MW freq sweep, we tell the PB to "WAIT" between each averaging block for each freq
        and tell it to resume after adjusting the MW freq

        The second way is significantly faster for a slow DAQ (e.g. USB DAQ). It shouldn't do much for a PCI DAQ which has much less read/write overhead
        :param in_data: I don't know what this is
        :return: Nothing
        """

        self.define_sweep_parameters()
        self._configure_param_array()
        self._configure_instruments_start_of_script()

        self.data['params'] = self.params
        self.data['value'] = []

        # Divides the total number of averages requested into a number of slices of MAX_AVERAGES_PER_SCAN and a remainder.
        # This is required because the pulseblaster won't accept more than ~4E6 loops (22 bits available to store loop
        # number) so need to break it up into smaller chunks (use 1E6 so initial results display faster)

        self.pulse_sequences, self.tau_list, self.measurement_gate_width = self.create_pulse_sequences()
        self.num_averages = self.settings['num_averages']

        # Keeps track of index of current pulse sequence for plotting
        self.sequence_index = 0

        for loop_index in range(self.num_averages):
            time_start = time.time()
            if breaker:
                break

            param_indices = list(range(len(self.params)))
            if self.settings['randomize']:
                np.random.shuffle(param_indices)

            for i in param_indices:
                if self._abort:
                    print('Aborting!!')
                    if self.instruments['PB']['instance'].settings['PB_type'] == 'USB':
                        self.instruments['PB']['instance'].stop_pulse_seq()
                    elif self.instruments['PB']['instance'].settings['PB_type'] == 'PCI':
                        try:
                            if self.instruments['PB']['instance'].pb.pb_read_status() & 0b100 == 0b100:  # Check if PB is running
                                self.instruments['PB']['instance'].pb.pb_close()
                            print('Sent close command to PB')
                        except:
                            print('Unable to read PB status. Did not send close command.')

                    self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                    self.log('Aborted pulseblaster script during loop')
                    breaker = True
                    break

                self.param_current = self.params[i]
                self._configure_instruments_start_of_sweep(self.param_current)
                self._configure_subscripts_start_of_sweep(self.param_current)
                time.sleep(self.sweep_params['param_switching_time'])

                self.scripts['backaction_pulse_seq'].run()
                self.scripts['lock_in_daq_read'].run()

                """
                Process daq data here
                """


            #  Save data after each averaging block so that we can analyze data externally on the go
            if self.settings['save']:
                self.save_data()

            time_elapsed = time.time() - time_start
            self.log("Completed loop %i of %i in %s" %
                     (loop_index + 1, self.num_averages, str(datetime.timedelta(seconds=time_elapsed))[:-7]))

            # This is needed to update the plotting. We put it here after performing the fitting
            self.updateProgress.emit(float(loop_index + 1) / int(self.num_1E5_avg_pb_programs) * 100)

        self._configure_instruments_end_of_script()
        self.updateProgress.emit(float(loop_index + 1) / int(self.num_1E5_avg_pb_programs) * 100)


class BackactionMwFreqSweep(BackactionGeneric):
    """
    Sweeps the MW freq of the pi-pulse train and measures the mechanical response using the lock-in DAQ
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('freq_start', 2.87e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 1e8, float, 'end frequency of scan in Hz'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 50, int, 'number of frequencies in scan in Hz'),
        Parameter('mw_generator_switching_time', .01, float, 'time wait after switching center frequencies on generator (s)')
    ]

    _INSTRUMENTS = {'mw_gen': MicrowaveGenerator}

    _SCRIPTS = {'pi_pulse_train_backaction': PiPulseTrainBackaction, 'lock_in_daq_read': LockInDaqRead}
