
# This file is part of pylabcontrol, software for laboratory equipment control for scientific experiments.
# Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell
#
# pylabcontrol is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pylabcontrol is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.


# this is the gui for the measurment pc

import ctypes
import sys

from PyQt5 import QtGui, QtWidgets

# from .pylabcontrol.gui.gui_future_b26_toolkit import b26_gui
from b26_toolkit.gui import b26_gui

#work around to change taskbar icon
#myappid = 'lukinlab.b26.pythonlab' # arbitrary string
#ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModeFlID(myappid)

def run_gui(fname = 'C:\\Users\Experiment\PycharmProjects\\user_data\pythonlab_config.b26'):
    app = QtWidgets.QApplication(sys.argv)
    # fname = 'C:\\Users\\NV Experiment\\PycharmProjects\\user_data\\new_config.b26'
    fname = ' '

    try:
        print("pylabcontrol  Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell \n\n \
pylabcontrol is free software: you can redistribute it and/or modify \n \
it under the terms of the GNU General Public License as published by \n \
the Free Software Foundation, either version 3 of the License, or \n \
(at your option) any later version. \n\n \
pylabcontrol is distributed in the hope that it will be useful, \n \
but WITHOUT ANY WARRANTY; without even the implied warranty of \n \
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the \n \
GNU General Public License for more details. \n\n \
You should have received a copy of the GNU General Public License \n \
along with pylabcontrol.  If not, see <http://www.gnu.org/licenses/>.")

        ex = b26_gui.ControlMainWindowB26(fname)

        app.setWindowIcon(QtGui.QIcon('magnet_and_nv.ico'))

        ex.show()
        ex.raise_()
        sys.exit(app.exec_())

    except ValueError as e:

        if not e.message in['No config file was provided. abort loading gui...', '']:
            raise e


if __name__ == '__main__':
    run_gui()