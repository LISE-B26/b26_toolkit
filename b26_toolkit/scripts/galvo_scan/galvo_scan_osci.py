from b26_toolkit.scripts.galvo_scan.galvo_scan_generic import GalvoScanGeneric
from pylabcontrol.core import Script, Parameter
import numpy as np
from b26_toolkit.scripts import SetLaser
from b26_toolkit.instruments import RigolOscilloscope
import time
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from b26_toolkit.plotting.plots_1d import plot_psd
from pylabcontrol.data_processing.signal_processing import power_spectral_density

class GalvoScanOsci(GalvoScanGeneric):
    _DEFAULT_SETTINGS = [
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate'),
                   Parameter('y', 0, float, 'y-coordinate')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate'),
                   Parameter('y', 1.0, float, 'y-coordinate')
                   ]),
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                               corner: pta and ptb are diagonal corners of rectangle.\n \
                                                               center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 64, int, 'number of x points to scan'),
                   Parameter('y', 64, int, 'number of y points to scan')
                   ]),
        Parameter('time_per_pt', .002, [.0005, .001, .002, .005, .01, .015, .02, .05, .08, .1],
                  'time in s to measure at each point'),
        Parameter('offset', 0.0, float, 'voltage offset [V]'),
        Parameter('vert_scale', 0.001, float, 'voltage scale [V]'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'],
                  'return to the corn'),
        Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    _INSTRUMENTS = {'oscope': RigolOscilloscope}

    _SCRIPTS = {'SetLaser': SetLaser}

    def __init__(self, name=None, settings=None, instruments=None, scripts=None, log_function=None, data_path=None):
        '''
        Initializes GalvoScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        Script.__init__(self, name=name, settings=settings, instruments=instruments, scripts=scripts, log_function=log_function, data_path=data_path)
        self.oscope = self.instruments['oscope']['instance']

    def setup_scan(self):
        """
        setup the scan, i.e. identify the instruments and set up sample rate and such
        """
        self.oscope.settings['vert_scale'] = self.settings['vert_scale']
        self.oscope.update(self.oscope.settings)
        self.oscope.settings['offset'] = self.settings['offset']
        self.oscope.update(self.oscope.settings)
        self.oscope.set_auto_trigger(self.settings['time_per_pt'] / 10.0)
        self.scripts['SetLaser'].settings['daq_type'] = 'PCI'
        self.data['voltage'] = []

    def read_line(self, y_pos):
        set_laser_script = self.scripts['SetLaser']
        set_laser_script.settings['point']['y'] = y_pos

        line_data = np.zeros(len(self.x_array))
        for i in range(len(self.x_array)):
            if self._abort:
                break

            set_laser_script.settings['point']['x'] = self.x_array[i]
            set_laser_script.run()
            time.sleep(self.settings['time_per_pt'])
            line_data[i] = self.oscope.get_std_voltage()

        return line_data

    def get_galvo_location(self):
        """
        Returns the current position of the galvo. Requires a daq with analog inputs internally routed to the analog
        outputs (ex. NI6259. Note that the cDAQ does not have this capability).
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        point = self.scripts['SetLaser'].settings['point']
        return [point['x'], point['y']]

    def set_galvo_location(self, galvo_position):
        """
        sets the current position of the galvo
        galvo_position: list with two floats, which give the x and y position of the galvo mirror
        """
        if galvo_position[0] > 1 or galvo_position[0] < -1 or galvo_position[1] > 1 or galvo_position[1] < -1:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 1 V')

        self.scripts['SetLaser'].settings['point']['x'] = galvo_position[0]
        self.scripts['SetLaser'].settings['point']['y'] = galvo_position[1]

# class GalvoScanOsci(GalvoScanGeneric):
#     """
#     GalvoScan uses the apd, daq, and galvo to sweep across voltages while counting photons at each voltage,
#     resulting in an image in the current field of view of the objective.
#     """
#
#     _DEFAULT_SETTINGS = [
#         Parameter('DAQ_channels',
#                    [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
#                     Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
#                     Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1'], 'Daq channel used for counter')
#                   ])
#     ]
#
#     _SCRIPTS = {'setlaser': SetLaser, 'get_trace': KeysightOsciGetTimeTrace}
#
#     _ACQ_TYPE = 'point'
#
#     def __init__(self, scripts, name=None, settings=None, log_function=None, data_path=None):
#         '''
#         Initializes GalvoScan script for use in gui
#
#         Args:
#             instruments: list of instrument objects
#             name: name to give to instantiated script object
#             settings: dictionary of new settings to pass in to override defaults
#             log_function: log function passed from the gui to direct log calls to the gui log
#             data_path: path to save data
#
#         '''
#         self._DEFAULT_SETTINGS = self._DEFAULT_SETTINGS + GalvoScanGeneric._DEFAULT_SETTINGS
#         Script.__init__(self, name, settings=settings, scripts=scripts, log_function=log_function, data_path = data_path)
#
#     def setup_scan(self):
#
#         self.clockAdjust = int((self.settings['time_per_pt'] + self.settings['settle_time']) / self.settings['settle_time'])
#
#         [xVmin, xVmax, yVmax, yVmin] = self.pts_to_extent(self.settings['point_a'],self.settings['point_b'],self.settings['RoI_mode'])
#
#
#         self.x_array = np.linspace(xVmin, xVmax, self.settings['num_points']['x'], endpoint=True)
#         self.y_array = np.linspace(yVmin, yVmax, self.settings['num_points']['y'], endpoint=True)
#
#
#     def get_galvo_location(self):
#         """
#         returns the current position of the galvo
#         Returns: list with two floats, which give the x and y position of the galvo mirror
#         """
#         galvo_position = self.scripts['setlaser'].get_galvo_position()
#         return galvo_position
#
#     def set_galvo_location(self, galvo_position):
#         """
#         sets the current position of the galvo
#         galvo_position: list with two floats, which give the x and y position of the galvo mirror
#         """
#
#         if isinstance(galvo_position, list):
#             galvo_position = {'x':galvo_position[0], 'y':galvo_position[1]}
#         elif isinstance(galvo_position, dict):
#             pass
#         else:
#             print(('asdasdad galvo_position', galvo_position))
#             raise TypeError
#
#         self.scripts['setlaser'].update({'point':galvo_position}) #update position for laser pointer
#         self.scripts['setlaser'].run() # point laser
#
#     def read_point(self, x_pos, y_pos):
#         """
#         reads a point of data from the Oscilloscope
#         Args:
#             x_pos: x position of the scan
#             y_pos: y position of the scan
#
#         Returns:
#
#         """
#
#         # point laser to new position
#         galvo_position = {'x':float(x_pos), 'y':float(y_pos)}
#         self.set_galvo_location(galvo_position)
#
#
#         # acquire timetrace with osci
#         self.scripts['get_trace'].run()
#
#         # JG: keep a record of the metadata so that we know the timestep later, this is not a very elegant but the best I came up with for now
#         self.data['meta_data'] = self.scripts['get_trace'].data['meta_data']
#
#         return self.scripts['get_trace'].data['voltage']
#
#
#     def _plot(self, axes_list, data=None):
#         """
#         Plots the galvo scan image
#         Args:
#             axes_list: list of axes objects on which to plot the galvo scan on the first axes object
#             data: data (dictionary that contains keys image_data, extent) if not provided use self.data
#         """
#
#         print(('axes list lenght plot', len(axes_list)))
#         if data is None:
#             data = self.data
#
#         plot_fluorescence_new(data['image_data'], data['extent'], axes_list[0])
#
#         last_data = self.data['point_data'][-1]
#         dt = self.data['meta_data']['xincrement']
#         freq, psd = power_spectral_density(last_data, dt)
#         plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')
#
#
#     def _update_plot(self, axes_list):
#         """
#         updates the galvo scan image
#         Args:
#             axes_list: list of axes objects on which to plot plots the esr on the first axes object
#         """
#         print(('axes list lenght update', len(axes_list)))
#
#         update_fluorescence(self.data['image_data'], axes_list[0])
#
#         last_data = self.data['point_data'][-1]
#         dt = self.data['meta_data']['xincrement']
#         freq, psd = power_spectral_density(last_data, dt)
#         plot_psd(freq, psd, axes_list[1], y_scaling='log', x_scaling='log')
#         axes_list[1].hold(False)
#
#
#     def get_axes_layout(self, figure_list):
#         """
#         returns the axes objects the script needs to plot its data
#         the default creates a single axes object on each figure
#         This can/should be overwritten in a child script if more axes objects are needed
#         Args:
#             figure_list: a list of figure objects
#         Returns:
#             axes_list: a list of axes objects
#
#         """
#
#         axes_list = []
#         if self._plot_refresh is True:
#             for fig in figure_list:
#                 fig.clf()
#                 axes_list.append(fig.add_subplot(111))
#
#         else:
#             for fig in figure_list:
#                 axes_list.append(fig.axes[0])
#
#         return axes_list