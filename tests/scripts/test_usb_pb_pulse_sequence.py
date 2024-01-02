from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_generic import PulsedExperimentGeneric
from b26_toolkit.instruments.pulse_blaster import Pulse

from pylabcontrol.core import Parameter


class TestUsbPbPulseSequence(PulsedExperimentGeneric):

    _DEFAULT_SETTINGS = [Parameter('num_averages', 100000, int, 'number of averages'),]

    def _create_pulse_sequences(self):
        tau = 50.
        measurement_time = 250.
        laser_off_time = 1000.
        nv_reset_time = 1750.
        delay_readout = 30.
        delay_mw_readout = 100.


        pulse_collection = [[Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
                             Pulse('apd_readout', laser_off_time + tau + 2*40 + delay_readout, measurement_time),
                             Pulse('microwave_q',
                                   laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time,
                                   tau),
                             Pulse('laser',
                                   laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout,
                                   nv_reset_time),
                             Pulse('apd_readout',
                                   laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
                                   measurement_time)
                             ]]

        print(pulse_collection)

        # pulse_collection = [[Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
        #                      Pulse('apd_readout', laser_off_time + tau + 2*40 + 20, measurement_time),
        #                      Pulse('microwave_q',
        #                            laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time,
        #                            tau),
        #                      Pulse('laser',
        #                            laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout,
        #                            nv_reset_time),
        #                      Pulse('apd_readout',
        #                            laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
        #                            measurement_time)
        #                      ]]

        tau_list = [tau]

        return pulse_collection, tau_list, measurement_time