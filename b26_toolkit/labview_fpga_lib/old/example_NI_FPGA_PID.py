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

"""
Created on Jan 29 2016

@author: Jan Gieseler

# example code for how to use PI loop on FPGA, loop terminates after a few loop iteractions
"""

import time

import matplotlib.pyplot as plt

from src import old_lib as NI

if __name__ == '__main__':
    # ============= DEFINE PARAMETERS ==================================
    # ==================================================================
    parameters_PI = {
        "sample_period_PI" :4e5,
        "gains" : {'proportional': 1.0, 'integral':0.1},
        "setpoint" : 0,
        "piezo" : 0
    }


    # ===========================================
    #  === helper old_functions
    # ===========================================
    def makeFig():
        # plt.scatter(xList,yList) # I think you meant this
        time = list(range(len(piezo_list)))
        axarr[0].scatter(time, piezo_list)
        axarr[1].scatter(time, detector_list)


    # ============= RUN PI LOOP ========================================
    # ==================================================================
    # connect to FPGA and start it
    fpga = NI.NI7845R()
    fpga.start()

    # create PI (proportional integral) controler object
    PI = NI.NI_FPGA_PI(fpga, **parameters_PI)


    # set piezo to mean value
    PI.piezo = 15000

    # turn feedback on
    PI.status_PI = False

    # create figure
    plt.ion() # enable interactivity
    fig=plt.figure() # make a figure
    f, axarr = plt.subplots(2, sharex=True)
    # axarr[0].set_title('time')
    # axarr[0].set_xlabel('time')
    axarr[0].set_xlabel('time')
    axarr[0].set_ylabel('piezo signal (bits)')
    axarr[1].set_ylabel('detector signal (bits)')

    piezo_list=list()
    detector_list=list()
    plt.hold(True)

    # plot signal
    for i in range(20):
        signal = PI.detector
        piezo = PI.piezo
        time.sleep(0.5)
        print(('detector value {:d}'.format(signal)))
        print(('piezo value {:d}'.format(piezo)))

        piezo_list.append(piezo)
        detector_list.append(signal)


        makeFig()
        plt.draw()
        # todo: find a better timer  function similar to labview "wait until next ms multiple"
        # to take into account the execution time of the resto of the code in the loop
        plt.pause(0.2)

    # ============= STOP FPGA ==========================================
    # ==================================================================
    PI.status_PI = False
    fpga.stop()

    print((fpga.session, fpga.status))
    eval(input("Please type enter to exit and close plot..."))