from pylabcontrol.gui.windows_and_widgets import LoadDialog

class LoadDialogB26(LoadDialog):
    """
This class builds on the Loaddialog and adds more options for script iterators
    """
    def __init__(self, elements_type, elements_old={}, filename=''):
        super(LoadDialogB26, self).__init__(elements_type, elements_old=elements_old, filename=filename)
        self.cmb_looping_variable.addItems(['Iter NVs', 'Iter Points', 'Iter test'])


if __name__ == '__main__':
    import sys
    from PyQt4 import QtGui

    app = QtGui.QApplication(sys.argv)
    # ex = LoadDialog(elements_type = 'instruments', elements_old=instuments, filename="Z:\Lab\Cantilever\Measurements\\__tmp\\test.b26")
    # ex = LoadDialog(elements_type='scripts', elements_old=instuments)
    ex = LoadDialogB26(elements_type='scripts', filename='/Users/rettentulla/Projects/Python/user_data')


    ex.show()
    ex.raise_()

    if ex.exec_():
        values = ex.getValues()

    sys.exit(app.exec_())




#
# if __name__ == '__main__':
#     from PyQt4 import QtGui
#     import sys
#
#     app = QtGui.QApplication(sys.argv)
#
#     dialog = LoadDialogB26(elements_type="scripts")
#
#     app.setWindowIcon(QtGui.QIcon('magnet_and_nv.ico'))
#
#     dialog.show()
#     dialog.raise_()
#     sys.exit(app.exec_())
#
#
#     print('-----', dialog.vars() )
#
#     if dialog.exec_():
#         xx = str(dialog.txt_probe_log_path.text())
#     scripts = dialog.getValues()
#
#     print('asdsadadas', xx)