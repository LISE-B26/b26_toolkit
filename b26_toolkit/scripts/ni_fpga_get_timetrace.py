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

from b26_toolkit.instruments import NI7845RReadFifo
from pylabcontrol.core import Script, Parameter


class LabviewFpgaTimetrace(Script):
    # COMMENT_ME

    _DEFAULT_SETTINGS = [
        Parameter('dt', 200, int, 'sample period of acquisition loop in ticks (40 MHz)'),
        Parameter('N', 10000, int, 'numer of samples'),
        # Parameter('TimeoutBuffer', 0, int, 'time after which buffer times out in clock ticks (40MHz)'),
        Parameter('BlockSize', 1000, int, 'block size of chunks that are read from FPGA'),
    ]

    _INSTRUMENTS = {'fpga' : NI7845RReadFifo}

    _SCRIPTS = {}

    def __init__(self, instruments, name = None, settings = None, log_function = None, data_path = None):
        Script.__init__(self, name, settings, instruments, log_function= log_function, data_path = data_path)

        self.data = deque()



    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        def calculate_progress(loop_index):
            progress = int(100.0 * loop_index / number_of_reads)
            return progress

        self._recording = True
        self.data.clear() # clear data queue


        # reset FIFO
        self.instruments['fpga'].stop_fifo()
        block_size = self.settings['BlockSize']
        # self.instruments['fpga'].update({'fifo_size' :block_size * 2})
        time.sleep(0.1)
        self.instruments['fpga'].start_fifo()
        time.sleep(0.1)
        number_of_reads = int(np.ceil(1.0 * self.settings['N'] / self.settings['BlockSize']))
        self.log('number_of_reads: {:d}'.format(number_of_reads))
        N_actual = number_of_reads * block_size
        if N_actual!=self.settings['N']:
            self.log('warning blocksize not comensurate with number of datapoints, set N = {:d}'.format(N_actual))
            self.settings.update({'N' : N_actual})
            time.sleep(0.1)

        # apply settings to instrument
        instr_settings = {
            'SamplePeriodsAcq' : self.settings['dt'],
            'ElementsToWrite' : self.settings['N']
            # 'TimeoutBuffer' : self.settings['TimeoutBuffer']
        }
        self.instruments['fpga'].update(instr_settings)
        time.sleep(0.1)
        self.instruments['fpga'].update({'Acquire' : True})

        ai1 = np.zeros(N_actual)
        ai2 = np.zeros(N_actual)
        i = 0
        while i < number_of_reads:
            elem_written = self.instruments['fpga'].ElementsWritten
            if elem_written >= block_size:
                data = self.instruments['fpga'].read_fifo(block_size)

                ai1[i * block_size:(i + 1) * block_size] = deepcopy(data['AI1'])
                ai2[i * block_size:(i + 1) * block_size] = deepcopy(data['AI2'])
                i += 1

                progress = calculate_progress(i)

                self.data.append({
                    'AI1' : ai1,
                    'AI2' : ai2,
                    'elements_remaining': data['elements_remaining']
                })

                self.updateProgress.emit(progress)

        self._recording = False

        # if self.settings['save']:
        #     self.save_b26()

    def _plot(self, axes_list, data = None):
        """
        Plots the timetrace taken with NI FPGA
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        axes = axes_list[0]

        if data is None:
            data = self.data

        if isinstance(data, deque):
            r = self.data[-1]['AI1']
        else:
            r = self.data['AI1']

        dt = self.settings['dt']/40e6

        time = dt * np.arange(len(r))
        if max(time)<1e-3:
            time *= 1e6
            xlabel = 'time (us)'
        elif max(time)<1e0:
            time *= 1e3
            xlabel = 'time (ms)'
        elif max(time)<1e3:
            xlabel = 'time (s)'
        elif max(time)<1e6:
            time *= 1e-3
            xlabel = 'time (ks)'
        axes.plot(time, r)
        axes.set_xlabel(xlabel)


if __name__ == '__main__':

    import time
    import numpy as np
    from copy import deepcopy

    fpga = NI7845RReadFifo()

    print((fpga.settings))

    # reset FIFO
    block_size = 2 ** 8

    N = 2 * block_size
    dt = 2000

    time.sleep(0.1)
    print('----stop-----')
    fpga.stop_fifo()
    print('----config-----')
    # fpga.update({'fifo_size': block_size * 2})
    print('----start-----')
    fpga.start_fifo()
    time.sleep(0.1)
    number_of_reads = int(np.ceil(1.0 * N / block_size))
    print(('number_of_reads', number_of_reads))
    N_actual = number_of_reads * block_size

    # apply settings to instrument
    instr_settings = {
        'SamplePeriodsAcq': dt,
        'ElementsToWrite': N
    }
    fpga.update(instr_settings)
    time.sleep(0.1)

    print('----------')
    print((fpga.settings))
    print('----------')

    print(('ElementsWritten: ', fpga.ElementsWritten))
    fpga.update({'Acquire': True})

    # time.sleep(1)
    print((fpga.settings))

    ai1 = np.zeros(N_actual)
    ai2 = np.zeros(N_actual)
    i = 0
    while i < number_of_reads:
        elem_written = fpga.ElementsWritten
        if elem_written >= block_size:
            data = fpga.read_fifo(block_size)
            # print(i, 'AI1', data['AI1'])
            print((i, 'elements_remaining', data['elements_remaining']))
            ai1[i * block_size:(i + 1) * block_size] = deepcopy(data['AI1'])
            ai2[i * block_size:(i + 1) * block_size] = deepcopy(data['AI2'])
            i += 1

        print(('-----', i, '------', 'elem_written', elem_written))

    print(ai1)
    print('------------------------------------------------')
    print(ai2)
