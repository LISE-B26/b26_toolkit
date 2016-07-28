from PyLabControl.src.core import Script
from b26_toolkit.src.scripts import ESR, Rabi


class ESRAndRabi(Script):
    """
    Does both an ESR experiment and a Rabi experiment on an NV, using the reference frequency from the esr data.
    """

    _DEFAULT_SETTINGS = []

    _INSTRUMENTS = {}

    _SCRIPTS = {'esr': ESR, 'rabi': Rabi}

    def __init__(self, scripts, name = None, settings = None, log_function = None, timeout = 1000000000, data_path = None):

        Script.__init__(self, name, scripts = scripts, settings=settings, log_function=log_function, data_path = data_path)


    def _function(self):

        self.scripts['esr'].run()

        if self.scripts['esr'].data['fit_params'] is not None:
            if len(self.scripts['esr'].data['fit_params']) == 4:
                self.rabi_frequency = self.scripts['esr'].data['fit_params'][2]
            elif len(self.scripts['esr'].data['fit_params']) == 6:
                self.rabi_frequency = self.scripts['esr'].data['fit_params'][4]
            else:
                raise RuntimeError('Could not get fit parameters from esr script')

            if self.rabi_frequency < self.scripts['esr'].settings['freq_start']:
                self.log('Resonance frequency found ({:0.2e}) was below esr sweep range, aborting rabi attempt'.format(self.rabi_frequency))
            elif self.rabi_frequency > self.scripts['esr'].settings['freq_stop']:
                self.log('Resonance frequency found ({:0.2e}) was above esr sweep range, aborting rabi attempt'.format(self.rabi_frequency))
            else:
                self.log('Starting rabi with frequency {:.4e} Hz'.format(self.rabi_frequency))
                self.scripts['rabi'].settings['mw_frequency'] = float(self.rabi_frequency)
                self.scripts['rabi'].run()
        else:
            self.log('No resonance frequency found skipping rabi attempt')


    def plot(self, axes_list):
        if self._current_subscript_stage['current_subscript'].name == 'esr':
            self.scripts['esr'].plot(axes_list)
        elif self._current_subscript_stage['current_subscript'].name == 'rabi':
            self.scripts['rabi'].plot(axes_list)