from PyLabControl.src.core import Script, Parameter
from src.scripts import StanfordResearch_ESR, Rabi


class ESR_AND_RABI(Script):
    """
    Does both an ESR experiment and a Rabi experiment on an NV, using the reference frequency from the esr data.
    """

    _DEFAULT_SETTINGS = []

    _INSTRUMENTS = {}

    _SCRIPTS = {'esr': StanfordResearch_ESR, 'rabi': Rabi}

    def __init__(self, scripts, name = None, settings = None, log_function = None, timeout = 1000000000, data_path = None):

        Script.__init__(self, name, scripts = scripts, settings=settings, log_function=log_function, data_path = data_path)


    def _function(self):

        self.scripts['esr'].run()

        if len(self.scripts['esr'].data['fit_params']) == 4:
            self.rabi_frequency = self.scripts['esr'].data['fit_params'][2]
        elif len(self.scripts['esr'].data['fit_params']) == 6:
            self.rabi_frequency = self.scripts['esr'].data['fit_params'][4]
        else:
            raise RuntimeError('Could not get fit parameters from esr script')

        if self.rabi_frequency > self.scripts['esr'].settings['freq_start']:
            self.lot('Resonance frequency found was below esr sweep range, aborting rabi attempt')
        elif self.rabi_frequency < self.scripts['esr'].settings['freq_stop']:
            self.log('Resonance frequency found was above esr sweep range, aborting rabi attempt')
        else:
            self.scripts['rabi'].settings['mw_frequency'] = self.rabi_frequency
            self.scripts['rabi'].run()


    def plot(self, axes_list):
        if self._current_subscript_stage['current_subscript'] == 'esr':
            self.scripts['esr'].plot(axes_list)
        elif self._current_subscript_stage['current_subscript'] == 'rabi':
            self.scripts['rabi'].plot(axes_list)