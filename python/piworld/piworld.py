#!/usr/bin/env python

from input import InputWidget
from project import ProjectWidget
from output import OutputWidget
import sys,os
from PyQt4 import QtGui,QtCore


class MainWindow(QtGui.QMainWindow):
  def __init__(self):
    QtGui.QMainWindow.__init__(self)
    icondir = os.path.dirname(__file__)+"/icons/"

    # Set up main window.
    self.resize(800, 600)
    self.setWindowTitle('piWorld')
    self.setWindowIcon(QtGui.QIcon(icondir+'pi-logo.png'))

    projectWidget = ProjectWidget(self)
    self.setCentralWidget(projectWidget)

    # Set up action objects.
    exitAction = QtGui.QAction(QtGui.QIcon(icondir+'silk/cancel.png'),
                               'Exit', self)
    exitAction.setShortcut('Ctrl+Q')
    exitAction.setStatusTip('Exit application')
    self.connect(exitAction, QtCore.SIGNAL('triggered()'),
                 QtCore.SLOT('close()'))

    newAction = QtGui.QAction(QtGui.QIcon(icondir+'silk/application_add.png'),
                                  'New', self)
    newAction.setShortcut('Ctrl+N')
    newAction.setStatusTip('New simulation')

    helpAction = QtGui.QAction(QtGui.QIcon(icondir+'silk/help.png'),
                                   'Help', self)
    helpAction.setShortcut('Ctrl+H')
    helpAction.setStatusTip('piDirector help')

    self.statusBar()

    # Set up menubar.
    menubar = self.menuBar()
    fileMenu = menubar.addMenu('&File')
    fileMenu.addAction(exitAction)
    fileMenu.addAction(newAction)
    editMenu = menubar.addMenu('&Edit')
    editMenu.addAction(helpAction)
    helpMenu = menubar.addMenu('&Help')
    helpMenu.addAction(helpAction)
    

    # Set up toolbar.
    toolbar = self.addToolBar('Exit')
    toolbar.addAction(newAction)
    toolbar.addAction(exitAction)
    toolbar.addSeparator()
    toolbar.addAction(helpAction)

app = QtGui.QApplication(sys.argv)
main = MainWindow()
main.show()
sys.exit(app.exec_())
