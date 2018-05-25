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

from pylabcontrol.core import Script, Parameter

class ScriptTest(Script):
    """
Minimal Example Script that has only a single parameter (execution time)
    """

    _DEFAULT_SETTINGS = [
        Parameter('execution_time', 0.1, float, 'execution time of script (s)'),
        Parameter('p1', 0.1, float, 'asihdad')
    ]

    _INSTRUMENTS = {}
    _SCRIPTS = {}

    def __init__(self, name=None, settings=None, log_function = None, data_path = None):
        """
        Example of a script
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings, log_function= log_function, data_path = data_path)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        import time
        time.sleep(self.settings['execution_time'])



if __name__ == '__main__':
    a = ScriptTest()

    print(a)