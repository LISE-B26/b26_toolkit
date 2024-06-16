from pylabcontrol.core import Parameter, Script
from b26_toolkit.scripts import StroboscopicReadout
from b26_toolkit.instruments import NI9263, WlmMonitorSiV

import numpy as np
import random
import time

class RedLaserSpectroscopy(Script):

    _DEFAULT_SETTINGS = [
        Parameter('sweep_param', 'voltage', ['voltage', 'frequency'], 'which parameter to sweep by'),
        Parameter('point_a', -3., float, 'start voltage/frequency of scan'),
        Parameter('point_b', 3., float, 'end voltage/frequency of scan'),
        Parameter('range_type', 'start_stop', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('num_points', 100, int, 'number of frequencies in scan'),
        Parameter('settle_time', 0.01, float, 'settle time between voltages [s]'),
        Parameter('sweep_type', 'up', ['up', 'down', 'randomize'], 'type of sweep'),
        Parameter('ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'analog output channel for setting laser wavelength')
    ]

    _INSTRUMENTS = {'daq_ao': NI9263, 'wlm': WlmMonitorSiV}
    _SCRIPTS = {'strobe_readout': StroboscopicReadout}

    _VOLT_MAX = 3.
    _VOLT_MIN = -3.

    _FREQ_MAX = 471.
    _FREQ_MIN = 470.

    _PLOTTING_DETUNING = 470.4

    def _get_sweep_array(self, max, min):
        '''

        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        point_a = self.settings['point_a']
        point_b = self.settings['point_b']
        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if point_a > point_b:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True



            sweep_values = np.linspace(point_a, point_b, self.settings['num_points'])

        elif self.settings['range_type'] == 'center_range':
            if point_a < 2 * point_b:
                self.log('end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True

                sweep_values = np.linspace(point_a - point_b / 2,
                                          point_a + point_b / 2, self.settings['num_points'])

        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        if np.min(sweep_values[0]) < min or np.max(sweep_values[-1]) > max:  # freq range of the SRS
            self.log('start or stop frequency out of bounds')
            self._abort = True

        return sweep_values

    def _function(self):
        sweep_param = self.settings['sweep_param']
        if sweep_param == 'voltage':
            sweep_array = self._get_sweep_array(self._VOLT_MAX, self._VOLT_MIN)
            self.data = {'counts': np.zeros(len(sweep_array)), 'voltage': sweep_array, 'freq': []}

        elif sweep_param == 'frequency':
            sweep_array = self._get_sweep_array(self._FREQ_MAX, self._FREQ_MIN)
            self.data = {'counts': np.zeros(len(sweep_array)), 'freq': sweep_array}
        else:
            self.log('incorrect sweep_param')
            return

        if self._abort:
            return

        daq_ao = self.instruments['daq_ao']['instance']
        wlm = self.instruments['wlm']['instance']
        readout = self.scripts['strobe_readout']

        ao_channel = self.settings['ao_channel']
        settle_time = self.settings['settle_time']

        indices = list(range(len(sweep_array)))
        if self.settings['sweep_type'] == 'randomize':
            random.shuffle(indices)
        elif self.settings['sweep_type'] == 'down':
            indices.reverse()

        try:
            if sweep_param == 'voltage':
                wlm.lock = False
                time.sleep(1)
                daq_ao.set_analog_voltages({ao_channel: sweep_array[indices[0]]})
            else:
                wlm.setpoint = sweep_array[indices[0]]
            time.sleep(1)
            for ind in indices:
                if self._abort:
                    break

                if sweep_param == 'voltage':
                    daq_ao.set_analog_voltages({ao_channel: sweep_array[ind]})
                else:
                    wlm.setpoint = sweep_array[ind]

                time.sleep(settle_time)

                if sweep_param == 'voltage':
                    self.data['freq'].append(wlm.frequency)
                readout.run()
                self.data['counts'][ind] = readout.data['counts'][0]

                self.updateProgress.emit(int(self.progress))

            if sweep_param == 'voltage':
                daq_ao.set_analog_voltages({ao_channel: 0})

        except RuntimeError:
            self.log('DAQ possible in use')

    def _plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if len(data['counts']) > 0:
            freqs = np.array(data['freq'])
            counts = np.array(data['counts']).flatten()
            axes_list[0].plot((freqs - self._PLOTTING_DETUNING) * 1e3, counts[:len(freqs)], linewidth=1.25)
            axes_list[0].set_xlabel('frequency [GHz]')
            axes_list[0].set_ylabel('[kCounts/s]')
            axes_list[0].set_title('NV spectroscopy (detuning from {} THz)'.format(self._PLOTTING_DETUNING))

            if self.settings['sweep_param'] == 'voltage':
                axes_list[1].plot((freqs - self._PLOTTING_DETUNING) * 1e3, data['voltage'][:len(freqs)], linewidth=1.25)
                axes_list[1].set_ylabel('voltage [V]')
                axes_list[1].set_xlabel('frequency [GHz]')

    def _update_plot(self, axes_list, data=None):
        if data is None:
            data = self.data

        if data:
            freqs = np.array(data['freq'])
            counts = np.array(data['counts']).flatten()
            axes_list[0].lines[0].set_xdata((freqs - self._PLOTTING_DETUNING) * 1e3)
            axes_list[0].lines[0].set_ydata(counts[:len(freqs)])
            axes_list[0].relim()
            axes_list[0].autoscale_view()

            if self.settings['sweep_param'] == 'voltage':
                axes_list[1].lines[0].set_ydata(data['voltage'][:len(freqs)])
                axes_list[1].lines[0].set_xdata((freqs - self._PLOTTING_DETUNING) * 1e3)
                axes_list[1].relim()
                axes_list[1].autoscale_view()