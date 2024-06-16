from pylabcontrol.core import Parameter, Instrument
from subprocess import Popen, PIPE, STDOUT
import time
import pyvisa as visa

class WlmMonitorSiV(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('server_name', 'wlm_monitor_server', str, 'server name'),
        Parameter('python_path', 'C:\\Users\\Experiment\\PycharmProjects\\pylabnet\\env\\Scripts\\python.exe', str, 'python path'),
        Parameter('launcher_path', 'C:\\Users\\Experiment\\PycharmProjects\\pylabnet\\pylabnet\\launchers\\launcher.py',str, 'launcher path'),
        Parameter('logip', '25.1.27.209', str, 'ip address of wlm monitor'),
        Parameter('logport', 12347, int, 'port of wlm monitor'),
        Parameter('script', 'wlm_monitor', str, 'pylabnet script name'),
        Parameter('num_clients', 32, int, 'number of clients'),
        Parameter('config', 'test', str, 'config filename'),
        Parameter('debug', 0, int, 'debug flag'),
        Parameter('server_debug', 0, int, 'server debug flag'),
        Parameter('lab_name', 'B26', str, 'lab name'),
        Parameter('initial_readlines', 14, int, 'intialy io readlines')
    ])

    def __init__(self, name=None, settings=None):

        super().__init__(name, settings)

        cmd = self.settings['python_path'] + ' '
        cmd += self.settings['launcher_path'] + ' '
        cmd += '--logip ' + self.settings['logip']
        cmd += ' --logport ' + str(self.settings['logport'])
        cmd += ' --script ' + self.settings['script']
        cmd += ' --num_clients ' + str(self.settings['num_clients'])
        cmd += ' --config ' + self.settings['config']
        cmd += ' --debug ' + str(self.settings['debug'])
        cmd += ' --server_debug ' + str(self.settings['server_debug'])
        cmd += ' --lab_name ' + self.settings['lab_name']
        self.sb = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT, shell=True, bufsize=0)
        for i in range(self.settings['initial_readlines']):
            self.sb.stdout.readline()

    def __del__(self):
        self.sb.kill()

    def _write(self, msg):
        self.sb.stdin.write('{}\r\n'.format(msg).encode())

    def _read(self):
        return self.sb.stdout.readline().decode().rstrip()


    @property
    def frequency(self):
        self._write('FREQUENCY?')
        return float(self._read())

    @property
    def setpoint(self):
        self._write('SETPOINT?')
        return float(self._read())


    @setpoint.setter
    def setpoint(self, sp):
        self._write('SETPOINT {}'.format(sp))
        if self._read() != 'ACK':
            raise IOError('no ACK received')

    @property
    def lock(self):
        self._write('LOCK?')
        return bool(int(self._read()))

    @lock.setter
    def lock(self, lock):
        self._write('LOCK {}'.format(int(lock)))
        if self._read() != 'ACK':
            raise IOError('no ACK received')

    @property
    def error_status(self):
        self._write('ERROR?')
        return bool(int(self._read()))

    def _PROBES(self):
        return {}

    def read_probes(self, key=None):
        return {}

class TLB6300LN(Instrument):
    _DEFAULT_SETTINGS = Parameter([
        Parameter('port', 'COM6', str, 'com port'),
        Parameter('baudrate', 19200, int, 'baud rate'),
        Parameter('setpoint', 637.1, float, 'wavelength setpoint controller in nm'),
        Parameter('wavelength', 637.1, float, 'sensed wavelength in nm')
    ])


    def __init__(self, name=None, settings=None):
        super().__init__(name, settings)

        rm = visa.ResourceManager()

        try:
            self.velocity = rm.open_resource('idk')
            self.velocity.baud_rate = self.settings['baudrate']
            self.velocity.query('*IDN?')
        except:
            self.log('unable to communicate with instrument')
            self._abort = True

    def update(self, settings):
        super().update(settings)

        if 'setpoint' in settings:
            setpoint = settings['setpoint']
        elif 'wavelength' in settings:
            setpoint = settings['wavelength']
        else:
            return

        self.velocity.write('SOURCE:WAVELENGTH {}'.format(str(setpoint)))


    @property
    def _PROBES(self):
        return {
            'wavelength': 'sensed wavelength in nm',
            'setpoint': 'wavelength setpoint controller in nm'
        }

    def read_probes(self, key=None):
        if key == 'setpoint':
            return self.velocity.query('SOURCE:WAVELENGTH?')
        elif key == 'wavelength':
            return self.velocity.query('SENSE:WAVELENGTH?')
        else:
            raise KeyError('unknown probe')

if __name__ == '__main__':
    wlm = WlmMonitorSiV()
    i = 0
    while True:
        print('get wavelength')
        print(wlm.frequency)
        time.sleep(5)
        i += 1

