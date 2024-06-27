"""
This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

pylabcontrol is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pylabcontrol is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.
"""


from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster
from pylabcontrol.core import Script, Parameter
#from b26_toolkit.scripts import FindNv
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric


class PulsedExperimentTracking(PulsedExperimentGeneric):
    """
    This class is a base class that should be inherited by all classes that utilize the pulseblaster for experiments. The
    _function part of this class takes care of high-level interaction with the pulseblaster for experiment control and optionally
    the daq for reading counter input (usually from the APD). It also provides all of the functionality needed to run a
    standard Script such as plotting.
    To use this class, the inheriting class need only overwrite _create_pulse_sequences to create the proper pulse sequence
    for a given experiment
    """

    _DEFAULT_SETTINGS = [
        Parameter('track_nv', [
            Parameter('below_threshold', False, bool, 'run FindNv whenever counts fall below a given threshold'),
            Parameter('threshold', 0.85, float, 'run FindNv whenever counts fall below this threshold'),
            Parameter('init_fluor', 20., float, 'initial fluorescence of the NV to compare to, in kcps'),
            Parameter('before_block', False, bool, 'run FindNv before each averaging block')]),
        Parameter('track_focus', [
            Parameter('threshold', 0.85, float, 'run FindNv whenever counts fall below this threshold'),
            Parameter('init_fluor', 20., float, 'initial fluorescence of the NV to compare to, in kcps'),
            Parameter('before_block', False, bool, 'run FindNv before each averaging block')])
    ]

    #_DEFAULT_SETTINGS = PulsedExperimentGeneric._DEFAULT_SETTINGS + _DEFAULT_SETTINGS
    _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster}

    # To enable tracking in a child script, simply enable the following script(s). The name given to the script (in parentheses, e.g. 'fin_nv') must be
    # the same as the following, but the script (e.g. FindNV) can be chosen as desired
    # Do NOT enable the tracking scripts here, enable them in the child scripts. Otherwise it will break some dependencies.
    #_SCRIPTS = {'find_nv': FindNV, 'autofocus': AutoFocusDaqMDT693A}
    #_SCRIPTS = {'find_nv': FindNv}

    def __init__(self, instruments, scripts, name=None, settings=None, log_function=None, data_path=None):
        """
        Standard script initialization
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        #self._DEFAULT_SETTINGS = self._DEFAULT_SETTINGS + PulsedExperimentGeneric._DEFAULT_SETTINGS
        self._DEFAULT_SETTINGS += PulsedExperimentGeneric._DEFAULT_SETTINGS
        self._DEFAULT_SETTINGS += PulsedExperimentTracking._DEFAULT_SETTINGS

        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.ref_index = 0

    def _track_nv(self, scenario, counts_temp=None):
        """
        Performs action to track to NV
        :param scenario:
                force -> run FindNv without checking any condition
                before_block -> run FindNv if 'before_block' in settings is True
                check_threshold -> run FindNv if counts are below threshold
        :param counts_temp: counts from the pulse sequence, to be checked against the threshold fluorescence, only needed if scenario is 'check_threshold'
        :return: whether to abort script, if FindNv fails consecutively multiple times
        """
        abort_script = False

        if scenario == 'before_block' and self.settings['track_nv']['before_block']:
            self.log('Running FindNv before averaging block')
            self.scripts['find_nv'].run(verbose=False)
            self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
        if scenario == 'force':
            self.log('Running FindNv on manual request')
            self.scripts['find_nv'].run(verbose=False)
            self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
        elif scenario == 'check_threshold' and self.settings['track_nv']['below_threshold']:
            threshold = self.settings['track_nv']['threshold']
            init_fluor = self.settings['track_nv']['init_fluor']
            findnv_attempts = 0
            counts_unsatisfactory = (1 + (1 - threshold)) * init_fluor < counts_temp or threshold * init_fluor > counts_temp
            self.log('Counts are below threshold; running FindNv')
            while counts_unsatisfactory:
                self.scripts['find_nv'].run(verbose=False)
                self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
                findnv_attempts += 1
                threshold *= 0.9
                counts_unsatisfactory = (1 + (1 - threshold)) * init_fluor < counts_temp or threshold * init_fluor > counts_temp
                if findnv_attempts >= 3:
                    self.log('Error: FindNv was unsuccessful after 5 attempts. Aborting script')
                    abort_script = True
                    break
                elif counts_unsatisfactory:
                    self.log('Counts still below threshold after FindNv, running FindNv again with lower threshold of %.1f kCt/s' % (threshold*init_fluor))

        return abort_script
