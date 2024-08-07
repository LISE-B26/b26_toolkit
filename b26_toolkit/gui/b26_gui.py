
from pylabcontrol.gui.windows_and_widgets.main_window import MainWindow
from b26_toolkit.gui.b26_load_dialog import LoadDialogB26
from PyQt5.QtWidgets import QLabel, QGraphicsOpacityEffect
from PyQt5.QtCore import pyqtSlot
import os
from pylabcontrol.core import Script
from playsound import playsound
from glob import glob
import random
from PyQt5 import QtGui

class ControlMainWindowB26(MainWindow):

    startup_msg = '\n\n\
    ======================================================\n\
    =============== Starting B26 toolkit  =============\n\
    ======================================================\n\n'

    def __init__(self, filename=None, argv=None):

        super(ControlMainWindowB26, self).__init__(filename)
        #self.setStyleSheet("background-color: #fae6fb;")
        #self.setStyleSheet("background-color: gray")
        # self.setStyleSheet("background-image: url(C:/Users/Experiment/Pictures/jojo.jpg); opacity: 0.01;")
        self.tree_scripts.setAlternatingRowColors(False)
        if argv is not None and len(argv) > 1:
            self.sound = argv[1] == 'True'
        else:
            self.sound = False

        self.setWindowIcon(QtGui.QIcon('poop.png'))


    def load_scripts(self, verbose=False):
        """
        opens file dialog to load scripts into gui
        """

        # update scripts so that current settings do not get lost
        for index in range(self.tree_scripts.topLevelItemCount()):
            script_item = self.tree_scripts.topLevelItem(index)
            self.update_script_from_item(script_item)

        dialog = LoadDialogB26(elements_type="scripts", elements_old=self.scripts,
                            filename=self.gui_settings['scripts_folder'])
        if dialog.exec_():
            self.gui_settings['scripts_folder'] = str(dialog.txt_probe_log_path.text())
            scripts = dialog.get_values()
            added_scripts = set(scripts.keys()) - set(self.scripts.keys())
            removed_scripts = set(self.scripts.keys()) - set(scripts.keys())

            if verbose:
                print(('load_scripts.scripts', scripts))
                print(('load_scripts.added_scripts', added_scripts))
                print(('load_scripts.removed_scripts', removed_scripts))

            if 'data_folder' in list(self.gui_settings.keys()) and os.path.exists(self.gui_settings['data_folder']):
                data_folder_name = self.gui_settings['data_folder']
            else:
                data_folder_name = None

            script_dict = {name: scripts[name] for name in added_scripts}

            if verbose:
                print(('load_scripts.script_dict', script_dict))
                print(('scripts, instruments', self.scripts, self.instruments))

            # create instances of new instruments/scripts
            self.scripts, loaded_failed, self.instruments = Script.load_and_append(
                script_dict=script_dict,
                scripts=self.scripts,
                instruments=self.instruments,
                log_function=self.log,
                data_path=data_folder_name)

            assert not 'self.scripts' == {}

            if verbose:
                print(('self.scripts', self.scripts))
            # delete instances of new instruments/scripts that have been deselected
            for name in removed_scripts:
                del self.scripts[name]

    @pyqtSlot()
    def script_finished(self):
        super(ControlMainWindowB26, self).script_finished()
        if self.sound:
            mp3s = glob('../media/sound/*.mp3')
            if len(mp3s) > 0:
                try:
                    playsound(random.choice(mp3s), block=False)
                except:
                    print('Could not play sound')
