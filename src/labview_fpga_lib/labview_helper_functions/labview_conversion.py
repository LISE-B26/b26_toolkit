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

some helper function for conversion of units and such


"""




def int_to_voltage(integer):
    return (10*integer)/32767.

def voltage_to_int(voltage):
    # TODO: make it work for arrays and lists
    return int((voltage * 32767)/10)

def time_to_buffersize(time, ticks=56):
    return int(time / (ticks*0.000000025))

def buffersize_to_time(size, ticks=56):
    return size * (ticks*0.000000025)

if __name__ == '__main__':
    print((voltage_to_int(2.4)))
