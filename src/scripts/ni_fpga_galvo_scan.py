"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
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

import datetime
from copy import deepcopy

import numpy as np

from b26_toolkit.src.instruments.labview_fpga import FPGA_GalvoScan
from src.plotting.plots_2d import  plot_fluorescence_new, update_fluorescence
from PyLabControl.src.core import Script, Parameter
from b26_toolkit.src.scripts.galvo_scan_generic import GalvoScanGeneric
from b26_toolkit.src.instruments.labview_fpga import volt_2_bit, NI7845RMain

class FPGA_GalvoScan(GalvoScanGeneric):
    """
    GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
    resulting in an image in the current field of view of the objective.
    """

    _DEFAULT_SETTINGS = [
    ]

    _INSTRUMENTS = {'NI7845RMain': NI7845RMain}

    def __init__(self, instruments, name=None, settings=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        self._DEFAULT_SETTINGS = self._DEFAULT_SETTINGS + GalvoScanGeneric._DEFAULT_SETTINGS
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function, data_path = data_path)

    def setup_scan(self):
        instr = self.instruments['NI7845RMain']['instance']
        instr_settings = self.instruments['NI7845RMain']['settings']['galvo_scan']

        [xVmin, xVmax, yVmax, yVmin] = volt_2_bit(self.data['extent'])

        Nx, Ny = self.settings['num_points']['x'], self.settings['num_points']['y']

        dVx = int((xVmax-xVmin) / Nx)
        dVy = int((yVmax - yVmin) / Ny)
        meas_per_pt = int(self.settings['time_per_pt'] * 400e3) # sample frequency is 400kHz
        settle_time = int(self.settings['settle_time'] * 1e6) # settle time of FPGA is in us

        instr_settings.update({
            'Vmin_x':xVmin,
            'Vmin_y':yVmin,
            'dVmin_x':dVx,
            'dVmin_y':dVy,
            'Nx':Nx,
            'Ny':Ny,
            'meas_per_pt': meas_per_pt,
            'settle_time':settle_time
        })


        instr.update({'galvo_scan':instr_settings})

        instr.stop_fifo()

        instr.start_fifo()

        started = instr.set_run_mode('galvo')

        if started is False:
            self.log('WARNING: galvo subvi did not start!!')

    def get_galvo_location(self):
        """
        returns the current position of the galvo
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        galvo_position = self.instruments['NI7845RMain']['instance'].AO0, self.instruments['NI7845RMain']['instance'].AO1
        return galvo_position

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        if galvo_position[0] > 1 or galvo_position[0] < -1 or galvo_position[1] > 1 or galvo_position[1] < -1:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 1 V')

        pt = volt_2_bit(galvo_position)

        self.instruments['NI7845RMain']['instance'].AO0 = pt[0]
        self.instruments['NI7845RMain']['instance'].AO1 = pt[1]


    def read_line(self, y_pos):
        """
        reads a line of data from the DAQ
        Args:
            y_pos: y position of the scan

        Returns:

        """
        instr = self.instruments['NI7845RMain']['instance']

        Nx = self.settings['num_points']['x']

        print('asdadsaada reading line ', Nx)
        line_data = instr.read_fifo(Nx)


        print('====== gggggg>>>>>>', line_data)
        return line_data['signal']

if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'FPGA_GalvoScan': 'FPGA_GalvoScan'})

    print('script',script)
    print('failed', failed)
    gs = script['FPGA_GalvoScan']
    print(gs)

    print(gs.instruments['FPGA_GalvoScan']['instance'].settings)

    gs.run()
    print(gs.data)

