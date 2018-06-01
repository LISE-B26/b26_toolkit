"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit  is free software: you can redistribute it and/or modify
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

from pylabcontrol.core import Instrument, Parameter
from b26_toolkit.labview_fpga_lib.labview_fpga_error_codes import LabviewFPGAException
import time
def volt_2_bit(volt):
    """
    converts a voltage value into a bit value
    Args:
        volt:

    Returns:

    """
    def convert(x):
        res = int(x / 10. * 32768.)
        res = min(res,32767) # cap at max value 32767
        res = max(res, -32768)  # cap at min value -32768

        return res

    if isinstance(volt, (int, float)):
        bit =  convert(volt)
    else:
        # convert to numpy array in case we received a list
        bit = [convert(x) for x in volt]
    return bit

def bit_2_volt(bit):
    """
    converts a voltage value into a bit value
    Args:
        bit:

    Returns:

    """
    if isinstance(bit, (int, float)):
        volt = bit * 10. / 32768
    else:
        # convert to numpy array in case we received a list
        volt = [int(x * 10./ 32768.) for x in bit]
    return volt

def seconds_to_ticks(seconds, clock_speed = 40e6):
    """
    convert seconds to ticks
    Args:
        seconds: time in seconds
        clock_speed: clock speed in Hz
    Returns:
        time in ticks
    """

    if isinstance(seconds, (int, float)):
        ticks = int(clock_speed * seconds)
    else:
        # convert to numpy array in case we received a list
        ticks = [int(clock_speed * x) for x in seconds]

    return ticks

# ==================================================================================
# simple fpga program that reads analog inputs and outputs
# ==================================================================================
class NI7845RMain(Instrument):

    import b26_toolkit.labview_fpga_lib.main.main as FPGAlib

    _DEFAULT_SETTINGS = Parameter([
        Parameter('read_io',[
            Parameter('AO0', 0.0, float, 'analog output channel 0 in volt'),
            Parameter('AO1', 0.0, float, 'analog output channel 1 in volt'),
            Parameter('AO2', 0.0, float, 'analog output channel 2 in volt'),
            Parameter('AO3', 0.0, float, 'analog output channel 3 in volt'),
            Parameter('AO4', 0.0, float, 'analog output channel 4 in volt'),
            Parameter('AO5', 0.0, float, 'analog output channel 5 in volt'),
            Parameter('AO6', 0.0, float, 'analog output channel 6 in volt'),
            Parameter('AO7', 0.0, float, 'analog output channel 7 in volt'),
            Parameter('DIO4', False, bool, 'digital output channel 4 on/off'),
            Parameter('DIO5', False, bool, 'digital output channel 5 on/off'),
            Parameter('DIO6', False, bool, 'digital output channel 6 on/off'),
            Parameter('DIO7', False, bool, 'digital output channel 7 on/off')
        ]),
        Parameter('galvo_scan',[
            Parameter('Vmin_x', 0, int, 'minimum voltage in x in bit'),
            Parameter('Vmin_y', 0, int, 'minimum voltage in y in bit'),
            Parameter('dVmin_x', 0, int, 'voltage step in x in bit'),
            Parameter('dVmin_y', 0, int, 'voltage step in y in bit'),
            Parameter('Nx', 0, int, 'number of voltage steps in x'),
            Parameter('Ny', 0, int, 'number of voltage steps in y'),
            Parameter('meas_per_pt', 1, int, 'number of measurements per point'),
            Parameter('settle_time',200, int,'wait time (us) between points to allow galvo to settle'),
            # Parameter('fifo_size', int(2 ** 12), int, 'size of fifo for data acquisition'),
            Parameter('scanmode_x', 'forward', ['forward', 'backward', 'forward-backward'],'scan mode (x) onedirectional or bidirectional'),
            Parameter('scanmode_y', 'forward', ['forward', 'backward'], 'direction of scan (y)'),
            Parameter('detector_mode', 'APD', ['APD', 'DC', 'RMS'], 'return mean (DC) or rms of detector signal')
        ]),
        Parameter('general', [
            Parameter('run_mode', 'idle', ['idle', 'galvo_scan', 'read_io'], 'select execution of subvi'),
            Parameter('count_ms', 1, int,'loop time of main loop, determines update frequency in idle and read_io mode')
        ])
    ])

    _PROBES = {
        'AI0': 'analog input channel 0 in bit',
        'AI1': 'analog input channel 1 in bit',
        'AI2': 'analog input channel 2 in bit',
        'AI3': 'analog input channel 3 in bit',
        'AI4': 'analog input channel 4 in bit',
        'AI5': 'analog input channel 5 in bit',
        'AI6': 'analog input channel 6 in bit',
        'AI7': 'analog input channel 7 in bit',
        'DIO0': 'digital input channel 0',
        'DIO1': 'digital input channel 1',
        'DIO2': 'digital input channel 2',
        'DIO3': 'digital input channel 3',
        'AO0': 'analog output channel 0 in bit',
        'AO1': 'analog output channel 1 in bit',
        'AO2': 'analog output channel 2 in bit',
        'AO3': 'analog output channel 3 in bit',
        'AO4': 'analog output channel 4 in bit',
        'AO5': 'analog output channel 5 in bit',
        'AO6': 'analog output channel 6 in bit',
        'AO7': 'analog output channel 7 in bit',
        'Nx': 'Nx',
        'Ny': 'Ny',
        # 'loop_time':'loop_time',
        # 'DMA_elem_to_write':'DMA_elem_to_write',
        'meas_per_pt': 'meas_per_pt',
        'settle_time': 'settle_time',
        'run_mode':'run_mode'
    }
    def __init__(self, name = None, settings = None):
        # start fpga
        self.fpga = self.FPGAlib.NI7845R()
        self.fpga.start()
        super(NI7845RMain, self).__init__(name, settings)
        # self.update(self.settings)

    def __del__(self):
        self.fpga.stop()

    def read_probes(self, key):
        if key is None:
            super(NI7845RMain, self).read_probes()
        else:
            assert key in list(self._PROBES.keys()), "key assertion failed %s" % str(key)
            value = getattr(self.FPGAlib, 'read_{:s}'.format(key))(self.fpga.session, self.fpga.status)
        return value

    def update(self, settings):
        super(NI7845RMain, self).update(settings)

        for key, value in settings.items():
            if key == 'read_io':
                for subkey, subvalue in value.items():
                    if subkey in ['AO0', 'AO1', 'AO2', 'AO3', 'AO4', 'AO5', 'AO6', 'AO7']:
                        getattr(self.FPGAlib, 'set_{:s}'.format(subkey))(volt_2_bit(subvalue), self.fpga.session, self.fpga.status)
                    elif subkey in ['DIO4', 'DIO5', 'DIO6', 'DIO7']:
                        getattr(self.FPGAlib, 'set_{:s}'.format(subkey))(subvalue, self.fpga.session, self.fpga.status)
                    else:
                        raise KeyError
            elif key == 'galvo_scan':
                for subkey, subvalue in value.items():
                    if subkey in ('scanmode_x', 'scanmode_y', 'detector_mode'):
                        subvalue = self._DEFAULT_SETTINGS.valid_values['galvo_scan'][subkey].index(subvalue)
                        getattr(self.FPGAlib, 'set_{:s}'.format(subkey))(subvalue, self.fpga.session, self.fpga.status)
                    elif subkey in ('Vmin_x', 'Vmin_y', 'dVmin_x', 'dVmin_y', 'Nx', 'Ny', 'meas_per_pt', 'settle_time'):
                        getattr(self.FPGAlib, 'set_{:s}'.format(subkey))(subvalue, self.fpga.session, self.fpga.status)
                    # elif subkey in ('fifo_size'):
                    #     pass
                    else:
                        raise KeyError
            elif key == 'general':
                for subkey, subvalue in value.items():
                    if subkey in ('run_mode'):
                        self.set_run_mode(subvalue)
                    elif subkey in ('count_ms'):
                        getattr(self.FPGAlib, 'set_{:s}'.format(subkey))(subvalue, self.fpga.session, self.fpga.status)
                    else:
                        raise KeyError

    def set_run_mode(self, mode, max_attempts=20):
        """
        start the the simple read io acquisition loop in the FPGA
        Returns:
            boolen that indecates wether start was successful or not
        """
        if isinstance(mode, str):
            if mode.lower() == 'read_io':
                mode = self.FPGAlib.READ_IO
            elif mode.lower() == 'idle':
                mode = self.FPGAlib.IDLE
            elif mode.lower() in ('galvo', 'galvo_scan', 'galvoscan'):
                mode = self.FPGAlib.GALVO_SCAN

        assert mode in (self.FPGAlib.READ_IO, self.FPGAlib.IDLE, self.FPGAlib.GALVO_SCAN)

        for i in range(max_attempts):

            getattr(self.FPGAlib, 'set_run_mode')(mode, self.fpga.session,self.fpga.status)

            # wait a little before checking of acquisition worked
            time.sleep(0.1)

            # run_mode = self.FPGAlib.read_run_mode(self.fpga.session,self.fpga.status)
            run_mode = self.run_mode
            print(('run_mode (', mode, ')', run_mode))
            print(('XXXX', run_mode, mode, run_mode == mode))
            started = run_mode == mode
            if started:
                # successfully started acquisition
                break
        if started == False:
            print(('starting FPGA (set mode to {:d}) failed after {:d} attempts!!!'.format(mode, max_attempts)))
            print(('current mode: {:d}'.format(self.read_probes('run_mode'))))


        return started

    def abort_acquire(self, max_attempts = 20):
        """
        stops the the simple read io acquisition loop in the FPGA and goes into the idle loop
        Returns:
            boolen that indecates wether start was successful or not
        """
        success = self.set_run_mode(self, self.FPGAlib.IDLE, max_attempts)
        # todo: JG: add abort button to FPGA and call here

        return success

    def start_fifo(self):
        self.FPGAlib.start_FIFO(self.fpga.session, self.fpga.status)

    def stop_fifo(self):
        self.FPGAlib.stop_FIFO(self.fpga.session, self.fpga.status)

    def read_fifo(self, block_size):
        '''
        read a block of data from the FIFO
        :return: data from channels AI1 and AI2 and the elements remaining in the FIFO
        '''

        print(('ssssssss sadsad block_size', block_size))
        fifo_data = self.FPGAlib.read_FIFO(block_size, self.fpga.session, self.fpga.status)
        if str(self.fpga.status.value) != '0':
            raise LabviewFPGAException(self.fpga.status)
        # print('fifo data', fifo_data)
        print('ss===>>>ssssss sadsad')
        if self.settings['galvo_scan']['detector_mode'] == 'APD':
            time_per_pt = self.settings['galvo_scan']['meas_per_pt']/400e3 # measurement rate is 400 kHz
            fifo_data['signal'] *=  int(1e3/time_per_pt) # convert to counts per second

        # todo: JG: scale to volts!
        # else:
        #     fifo_data['signal'] = np.array(bit_2_volt(fifo_data['signal']))
        return fifo_data
# ==================================================================================
# simple fpga program that reads analog inputs and outputs
# ==================================================================================
class NI7845RReadWrite(Instrument):

    # import pylabcontrol.labview_fpga_lib.read_ai_ao.read_ai_ao as FPGAlib
    import b26_toolkit.labview_fpga_lib.main.main as FPGAlib

    _DEFAULT_SETTINGS = Parameter([
        Parameter('AO0', 0.0, float, 'analog output channel 0 in volt'),
        Parameter('AO1', 0.0, float, 'analog output channel 1 in volt'),
        Parameter('AO2', 0.0, float, 'analog output channel 2 in volt'),
        Parameter('AO3', 0.0, float, 'analog output channel 3 in volt'),
        Parameter('AO4', 0.0, float, 'analog output channel 4 in volt'),
        Parameter('AO5', 0.0, float, 'analog output channel 5 in volt'),
        Parameter('AO6', 0.0, float, 'analog output channel 6 in volt'),
        Parameter('AO7', 0.0, float, 'analog output channel 7 in volt'),
        Parameter('DIO4', False, bool, 'digital output channel 4 on/off'),
        Parameter('DIO5', False, bool, 'digital output channel 5 on/off'),
        Parameter('DIO6', False, bool, 'digital output channel 6 on/off'),
        Parameter('DIO7', False, bool, 'digital output channel 7 on/off')
    ])

    _PROBES = {
        'AI0': 'analog input channel 0 in bit',
        'AI1': 'analog input channel 1 in bit',
        'AI2': 'analog input channel 2 in bit',
        'AI3': 'analog input channel 3 in bit',
        'AI4': 'analog input channel 4 in bit',
        'AI5': 'analog input channel 5 in bit',
        'AI6': 'analog input channel 6 in bit',
        'AI7': 'analog input channel 7 in bit',
        'DIO0': 'digital input channel 0',
        'DIO1': 'digital input channel 1',
        'DIO2': 'digital input channel 2',
        'DIO3': 'digital input channel 3',
        'AO0': 'analog output channel 0 in bit',
        'AO1': 'analog output channel 1 in bit',
        'AO2': 'analog output channel 2 in bit',
        'AO3': 'analog output channel 3 in bit',
        'AO4': 'analog output channel 4 in bit',
        'AO5': 'analog output channel 5 in bit',
        'AO6': 'analog output channel 6 in bit',
        'AO7': 'analog output channel 7 in bit'
    }
    def __init__(self, name = None, settings = None):
        super(NI7845RReadWrite, self).__init__(name, settings)

        # start fpga
        self.fpga = self.FPGAlib.NI7845R()
        self.fpga.start()
        self.update(self.settings)

    def __del__(self):
        self.fpga.stop()

    def read_probes(self, key):
        assert key in list(self._PROBES.keys()), "key assertion failed %s" % str(key)
        value = getattr(self.FPGAlib, 'read_{:s}'.format(key))(self.fpga.session, self.fpga.status)
        return value


    def update(self, settings):
        super(NI7845RReadWrite, self).update(settings)

        for key, value in settings.items():
            if key in ['AO0', 'AO1', 'AO2', 'AO3', 'AO4', 'AO5', 'AO6', 'AO7']:
                getattr(self.FPGAlib, 'set_{:s}'.format(key))(volt_2_bit(value), self.fpga.session, self.fpga.status)
            elif key in ['DIO4', 'DIO5', 'DIO6', 'DIO7']:
                getattr(self.FPGAlib, 'set_{:s}'.format(key))(value, self.fpga.session, self.fpga.status)

    def start_acquire(self, max_attempts=20):
        """
        start the the simple read io acquisition loop in the FPGA
        Returns:
            boolen that indecates wether start was successful or not
        """

        # self.stop_fifo()
        # # start fifo
        # self.start_fifo()

        # start scan
        for i in range(max_attempts):

            getattr(self.FPGAlib, 'set_run_mode')(self.FPGAlib.READ_IO, self.fpga.session,
                                              self.fpga.status)  #

            # wait a little before checking of acquisition worked
            time.sleep(0.1)

            run_mode = self.run_mode
            print(('run_mode (', self.FPGAlib.READ_IO, ')', run_mode))
            started = run_mode == self.FPGAlib.READ_IO
            if started:
                # successfully started acquisition
                break
        if started == False:
            print(('starting FPGA failed after {:d} attempts!!!'.format(max_attempts)))
            print((self.read_probes()))

        return True

    def stop_acquire(self, max_attempts = 20):
        """
        stops the the simple read io acquisition loop in the FPGA and goes into the idle loop
        Returns:
            boolen that indecates wether start was successful or not
        """
        getattr(self.FPGAlib, 'set_run_mode')(self.FPGAlib.IDLE, self.fpga.session, self.fpga.status) #

        return True

# ==================================================================================
# fpga program that performs a galvo scan
# ==================================================================================

class FPGA_GalvoScan(Instrument):

    # import b26_toolkit.pylabcontrol.labview_fpga_lib.galvo_scan.galvo_scan as FPGAlib
    import b26_toolkit.labview_fpga_lib.main.main as FPGAlib

    _DEFAULT_SETTINGS = Parameter([
        Parameter('point_a',
                  [Parameter('x', -0.4, float, 'x-coordinate'),
                   Parameter('y', -0.4, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 0.4, float, 'x-coordinate'),
                   Parameter('y', 0.4, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'corner', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                       corner: pta and ptb are diagonal corners of rectangle.\n \
                                                       center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 120, int, 'number of x points to scan'),
                   Parameter('y', 120, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', 0.25, [ 0.25, 0.5, 1.0, 2.0, 5.0, 10.0], 'time (ms) to measure at each point'),
        Parameter('settle_time', 0.1, [ 0.1, 0.2, 0.5, 1.0, 2.0], 'wait time (ms) between points to allow galvo to settle'),
        Parameter('fifo_size', int(2**12), int, 'size of fifo for data acquisition'),
        Parameter('scanmode_x', 'forward', ['forward', 'backward', 'forward-backward'], 'scan mode (x) onedirectional or bidirectional'),
        Parameter('scanmode_y', 'forward', ['forward', 'backward'], 'direction of scan (y)'),
        Parameter('detector_mode', 'APD', ['APD', 'DC', 'RMS'], 'return mean (DC) or rms of detector signal'),
        # Parameter('piezo', 0.0, float, 'output voltage for piezo'),
        # Parameter('piezo_gain', 10.0, [1.0, 7.5, 10.0, 15.0], 'gain factor for the piezo controller '),
        # Parameter('galvo_x', 0.0, float, 'output voltage for galvo channel x'),
        # Parameter('galvo_y', 0.0, float, 'output voltage for galvo channel y'),
    ])

    _PROBES = {
        # 'Detector Signal': 'detector signal (AI4)',
        'elements_written_to_dma' : 'elements written to DMA',
        # 'DMATimeOut': 'DMATimeOut',
        # # 'ix': 'ix',
        # 'iy':'iy',
        # # 'detector_signal':'detector_signal',
        # 'acquire':'acquire',
        'executing_subvi':'executing_subvi',
        'Nx': 'Nx',
        'Ny': 'Ny',
        # 'loop_time':'loop_time',
        # 'DMA_elem_to_write':'DMA_elem_to_write',
        'meas_per_pt':'meas_per_pt',
        'settle_time':'settle_time',
        # 'failed':'failed',
        # 'piezo':'piezo'
    }
    def __init__(self, name = None, settings = None):
        super(FPGA_GalvoScan, self).__init__(name, settings)

        # start fpga
        self.fpga = self.FPGAlib.NI7845R()
        self.fpga.start()
        self.update(self.settings)
        print('FPGA_GalvoScan initialized')

    def __del__(self):
        print(('stopping fpga {:s}'.format(self.name)))
        #define NiFpga_GalvoScan_Bitfile "C:\\Users\\Experiment\\PycharmProjects\\PythonLab\\pylabcontrol\\labview_fpga_lib\\galvo_scan\\NiFpga_GalvoScan.lvbitx" stop()

    def read_probes(self, key = None):

        if key is None:
            super(FPGA_GalvoScan, self).read_probes()
        else:
            assert key in list(self._PROBES.keys()), "key assertion failed %s" % str(key)
            # if key == 'ElementsWritten':
            #     key = 'elements_written_to_dma'
            value = getattr(self.FPGAlib, 'read_{:s}'.format(key))(self.fpga.session, self.fpga.status)
        return value

    @staticmethod
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

    # explicitely write piezo as property. I think this is not needed. try to remove later..
    @property
    def piezo(self):
        return bit_2_volt(self.read_probes('piezo'))
    @piezo.setter
    def piezo(self, value):
        """

        Args:
            value: piezo value in V this is the value output from the FPGA and will be amplified by a piezo controller!

        Returns:

        """
        self.update({'piezo':value})

    def update(self, settings):
        super(FPGA_GalvoScan, self).update(settings)

        for key, value in settings.items():
            if key in ['point_a', 'point_b', 'RoI_mode', 'num_points']:
                [xVmin, xVmax, yVmax, yVmin] = self.pts_to_extent(self.settings['point_a'],
                                                                  self.settings['point_b'],
                                                                  self.settings['RoI_mode'])
                # convert volt to I16 values
                [xVmin, xVmax, yVmax, yVmin] = volt_2_bit([xVmin, xVmax, yVmax, yVmin])

                Nx, Ny = self.settings['num_points']['x'], self.settings['num_points']['y']

                dVmin_x = int((xVmax- xVmin) / Nx)
                dVmin_y = int((yVmax - yVmin) / Ny)

                getattr(self.FPGAlib, 'set_Vmin_x')(xVmin, self.fpga.session, self.fpga.status)
                getattr(self.FPGAlib, 'set_Vmin_y')(yVmin, self.fpga.session, self.fpga.status)
                getattr(self.FPGAlib, 'set_Nx')(Nx, self.fpga.session, self.fpga.status)
                getattr(self.FPGAlib, 'set_Ny')(Ny, self.fpga.session, self.fpga.status)
                getattr(self.FPGAlib, 'set_dVmin_x')(dVmin_x, self.fpga.session, self.fpga.status)
                getattr(self.FPGAlib, 'set_dVmin_y')(dVmin_y, self.fpga.session, self.fpga.status)
            # elif key in ['piezo', 'piezo_gain']:
            #     v = volt_2_bit(self.settings['piezo']/self.settings['piezo_gain'])
            #     getattr(self.FPGAlib, 'set_piezo_voltage')(v, self.fpga.session, self.fpga.status)
            # elif key in ['galvo_x', 'galvo_x']:
            #     v = volt_2_bit(self.settings[key])
            #     getattr(self.FPGAlib, 'set_{:s}'.format(key))(v, self.fpga.session, self.fpga.status)
            elif key in ['scanmode_x', 'scanmode_y', 'detector_mode']:
                index = [i for i, x in enumerate(self._DEFAULT_SETTINGS.valid_values[key]) if x == value][0]
                getattr(self.FPGAlib, 'set_{:s}'.format(key))(index, self.fpga.session, self.fpga.status)
            elif key in ['settle_time']:
                settle_time = int(value*1e3)
                getattr(self.FPGAlib, 'set_settle_time')(settle_time, self.fpga.session, self.fpga.status)
            elif key in ['time_per_pt']:
                measurements_per_pt = int(value/0.25)
                getattr(self.FPGAlib, 'set_meas_per_pt')(measurements_per_pt, self.fpga.session, self.fpga.status)
            elif key in ['fifo_size']:
                actual_fifo_size = self.FPGAlib.configure_FIFO(value, self.fpga.session, self.fpga.status)
                # print('requested ', value )
                # print('actual_fifo_size ', actual_fifo_size)

    def start_fifo(self):
        self.FPGAlib.start_FIFO(self.fpga.session, self.fpga.status)

    def stop_fifo(self):
        self.FPGAlib.stop_FIFO(self.fpga.session, self.fpga.status)


    def start_acquire(self, max_attempts = 20):
        """
        started the acquisition loop in the FPGA
        max_attempts: maximum number of attempts to start acquisition
        Returns:
            boolen that indecates wether start was successful or not
        """

        self.stop_fifo()
        # start fifo
        self.start_fifo()


        # start scan
        for i in range(max_attempts):
            # getattr(self.FPGAlib, 'set_acquire')(True, self.fpga.session, self.fpga.status)
            getattr(self.FPGAlib, 'set_run_mode')(self.FPGAlib.GALVO_SCAN, self.fpga.session, self.fpga.status) #

            # wait a little before checking of acquisition worked
            time.sleep(0.1)
            started = self.executing_subvi
            if started:
                # successfully started acquisition
                break
        if started == False:
            print(('starting FPGA failed after {:d} attempts!!!'.format(max_attempts)))
            print((self.read_probes()))

        return started

    def abort_acquire(self):
        getattr(self.FPGAlib, 'set_stop_all')(True, self.fpga.session, self.fpga.status)

    def read_fifo(self, block_size):
        '''
        read a block of data from the FIFO
        :return: data from channels AI1 and AI2 and the elements remaining in the FIFO
        '''
        fifo_data = self.FPGAlib.read_FIFO(block_size, self.fpga.session, self.fpga.status)
        if str(self.fpga.status.value) != '0':
            raise LabviewFPGAException(self.fpga.status)
        # print('fifo data', fifo_data)

        if self.settings['detector_mode'] == 'APD':
            fifo_data['signal'] *=  int(1e3/self.settings['time_per_pt']) # convert to counts per second

        # todo: JG: scale to volts!
        # else:
        #     fifo_data['signal'] = np.array(bit_2_volt(fifo_data['signal']))
        return fifo_data



if __name__ == '__main__':
    import time
    from copy import deepcopy

    fpga = NI7845RMain()



    fpga.fpga.start()

    print((fpga.FPGAlib.setter_functions))
    print((fpga.settings))
    fpga.fpga.stop()





