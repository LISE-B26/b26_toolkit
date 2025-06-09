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
        Parameter('logport', 12352, int, 'port of wlm monitor'),
        Parameter('script', 'wlm_monitor', str, 'pylabnet script name'),
        Parameter('num_clients', 19, int, 'number of clients'),
        Parameter('config', 'test', str, 'config filename'),
        Parameter('debug', 0, int, 'debug flag'),
        Parameter('server_debug', 0, int, 'server debug flag'),
        Parameter('lab_name', 'B26', str, 'lab name'),
        Parameter('initial_readlines', 18, int, 'intialy io readlines')
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
    _SLEW_MIN = 0.01

    _DEFAULT_SETTINGS = Parameter([
        Parameter('visa_address', 'ASRL19::INSTR', str, 'com port'),
        Parameter('baudrate', 9600, int, 'baud rate'),
        Parameter('wavelength', 637.1, float, 'sensed wavelength in nm'),
        Parameter('start_wavelength', 636., float, 'starting wavelength for scan in nm'),
        Parameter('stop_wavelength', 638., float, 'stopping wavelength for scan in nm'),
        Parameter('fwd_scan_speed', 0.01, float, 'scanning speed in nm/s'),
        Parameter('rev_scan_speed', 0.01, float, 'scanning speed in nm/s')
    ])

    _COMMANDS = {
        'start_wavelength': ':WAVE:STAR ',
        'stop_wavelength': ':WAVE:STOP ',
        'fwd_scan_speed': ':WAVE:SLEW:FORW ',
        'rev_scan_speed': ':WAVE:SLEW:RET '
    }

    def __init__(self, name='None', settings=None):
        super().__init__(name, settings)
        self._is_connected = False

        rm = visa.ResourceManager()

        try:
            self._velocity = rm.open_resource(self.settings['visa_address'])
            self._is_connected = True
            self._velocity.write_termination = '\r'
            self._velocity.read_termination = '\r'
            self._velocity.baud_rate = self.settings['baudrate']

            self._WAVE_MAX = float(self._velocity.query(':WAVE? MAX')) - 0.1
            self._WAVE_MIN = float(self._velocity.query(':WAVE? MIN')) + 0.1
            self._SLEW_MAX = float(self._velocity.query(':WAVE:SLEW:MAX?')) - self._SLEW_MIN

            if self.settings['start_wavelength'] > self._WAVE_MAX \
                or self.settings['start_wavelength'] < self._WAVE_MIN \
                or self.settings['stop_wavelength'] > self._WAVE_MAX \
                or self.settings['stop_wavelength'] < self._WAVE_MIN:

                raise ValueError('start/stop wavelength value illegal')

            if self.settings['fwd_scan_speed'] < self._SLEW_MIN \
                or self.settings['fwd_scan_speed'] > self._SLEW_MAX \
                or self.settings['rev_scan_speed'] < self._SLEW_MIN \
                or self.settings['rev_scan_speed'] > self._SLEW_MAX:

                raise ValueError('slew rate value illegal')

            self.update(self.settings)
        except:
            self.log('unable to communicate with instrument')
            self._abort = True

    def __del__(self):
        if self.is_connected:
            self._velocity.close()

    def update(self, settings):
        super().update(settings)
        print(settings)
        for key, value in settings.items():
            if key in self._COMMANDS:
                self._velocity.write(self._COMMANDS[key] + str(value))
                msg = self._velocity.read()
                if msg != 'OK':
                    raise IOError('ack not received: {}'.format(msg))
            elif key == 'wavelength':
                self._velocity.write(':WAVE ' + str(value))
                msg = self._velocity.read()
                if msg != 'OK':
                    raise IOError('ack not received: {}'.format(msg))

                time.sleep(1)

                self._velocity.write(':OUTP:TRAC OFF')
                msg = self._velocity.read()
                if msg != 'OK':
                    raise IOError('ack not received: {}'.format(msg))

        ready = False
        while not ready:
            ready = bool(int(self._velocity.query('*OPC?')))

    @property
    def _PROBES(self):
        return {
            'wavelength': 'sensed wavelength in nm',
            'start_wavelength': 'starting wavelength for a scan in nm',
            'stop_wavelength': 'stopping wavelength for a scan in  nm',
            'fwd_scan_speed': 'forward scanning speed in nm/s',
            'rev_scan_speed': 'reverse scanning speed in nm/s'
        }

    def read_probes(self, key):
        assert (self._settings_initialized)  # will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())
        if key == 'wavelength':
            return float(self._velocity.query(':SENS:WAVE'))
        elif key in self._COMMANDS:
            return float(self._velocity.query(self._COMMANDS[key][:-1] + '?'))
        else:
            raise KeyError('unknown probe: {}'.format(key))

    def scan_finished(self):
        return bool(int(self._velocity.query('*OPC?')))

    def start_scan_wavelength(self):
        self._velocity.write(':OUTPUT:SCAN:RESET')
        msg = self._velocity.read()
        if msg != 'OK':
            raise IOError('ack not received: {}'.format(msg))

        self._velocity.write(':WAVE ' + str(self.settings['start_wavelength']))
        msg = self._velocity.read()
        if msg != 'OK':
            raise IOError('ack not received: {}'.format(msg))
        time.sleep(1)

        self._velocity.write(':OUTPUT:SCAN:START')
        msg = self._velocity.read()
        if msg != 'OK':
            raise IOError('ack not received: {}'.format(msg))

    def stop_scan_wavelength(self):
        self._velocity.write(':OUTPUT:SCAN:STOP')
        msg = self._velocity.read()
        if msg != 'OK':
            raise IOError('ack not received: {}'.format(msg))

    def set_full_range(self):
        self.update({'start_wavelength': self._WAVE_MIN,
                     'stop_wavelength': self._WAVE_MAX})

if __name__ == '__main__':
    # wlm = WlmMonitorSiV()
    # i = 0
    # while True:
    #     print('get wavelength')
    #     print(wlm.frequency)
    #     time.sleep(5)
    #     i += 1
    print('create red laser instance')
    redLaser = TLB6300LN()
    redLaser.start_scan_wavelength()
    while True:
        time.sleep(1)
        print(redLaser.wavelength)
    # redLaser.scan_wavelength()