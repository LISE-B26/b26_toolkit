from pylabcontrol.core import Script, Parameter
from b26_toolkit.instruments import B26PulseBlaster
from b26_toolkit.scripts import Daq_TimeTrace_NI6259

import time
import numpy as np

class Ringdown(Daq_TimeTrace_NI6259):
    _DEFAULT_SETTINGS = [
        Parameter('holdon', 1., float, 'time to hold on in s'),
        Parameter('holdoff', 1., float, 'time to hold off in s'),
        Parameter('pb_channel', 'atto_trig', str, 'channel to trigger drive'),
    ]

    _INSTRUMENTS = {'PB': B26PulseBlaster}
    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        self._DEFAULT_SETTINGS += super()._DEFAULT_SETTINGS
        self._INSTRUMENTS.update(super()._INSTRUMENTS)

        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        self.data = {'counts': []}

        sample_rate = float(1) / self.settings['integration_time']
        self.setup_daq(sample_rate)


        # maximum number of samples if total_int_time > 0
        if self.settings['acquisition_time'] > 0:
            number_of_samples = int(np.floor(self.settings['acquisition_time'] / self.settings['integration_time']))
        else:
            self.log('total measurement time must be positive. Abort script')
            return

        # self.progress = 50
        # self.updateProgress.emit(self.progress)
        self.data['start_time'] = time.time() # in seconds

        for i in range(self.settings['num_averages']):
            if self._abort:
                break

            self.log('starting run {} of {}'.format(i + 1, self.settings['num_averages']))
            daq_tasks = self.setup_daq_tasks(number_of_samples)

            # DS20191016 Record start time
            # print('Start time: {}'.format(time.time()))
            self.instruments['PB']['instance'].update({self.settings['pb_channel']: {'status': True}})
            time.sleep(self.settings['holdon'])
            if self._abort:
                break

            for daq, task in daq_tasks:
                # start task
                if daq is not None:
                    daq.run(task)
            time.sleep(self.settings['holdoff'])
            self.instruments['PB']['instance'].update({self.settings['pb_channel']: {'status': False}})


            data = self.read_daq_data(daq_tasks, sample_rate)

            self.data['counts'].append(data['counts'])


            for daq, task in daq_tasks:
                # start task
                if daq is not None:
                    daq.stop(task)

        self.data['counts'] = np.array(self.data['counts']).transpose()
