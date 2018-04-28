"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

"""
Created on Jan 29 2016

@author: Jan Gieseler

# example code for read and write analog inputs of labview FPGA
# connect AO0 (piezo out) to AI1 (detector in)
# Then run this example.
"""

import time

from src import old_lib as NI

if __name__ == '__main__':
    fpga = NI.NI7845R()

    print((fpga.session, fpga.status))
    fpga.start()

    print((fpga.session, fpga.status))



    AI = NI.AnalogInput(1,fpga)
    AO = NI.AnalogOutput(0,fpga)

    for i in range(0,20000,1000):
        AO.write(i)
        time.sleep(0.1)
        x = AI.read()
        print(('set {:05d}\t actual {:05d}\t error {:0.2f}%'.format(i, x, 100.* (x-i) / (x+i))))
    fpga.stop()

    print((fpga.session, fpga.status))

