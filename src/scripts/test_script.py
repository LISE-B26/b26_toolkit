from PyLabControl.src.core import Script, Parameter

class ScriptTest(Script):
    """
Minimal Example Script that has only a single parameter (execution time)
    """

    _DEFAULT_SETTINGS = [
        Parameter('execution_time', 0.1, float, 'execution time of script (s)')
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