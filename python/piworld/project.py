import sys,os,glob
from PyQt4 import QtGui,QtCore
from input import InputWidget
from output import OutputWidget

class ProjectWidget(QtGui.QWidget):
  def __init__(self, parent=None):
    QtGui.QWidget.__init__(self,parent)
    self.autoOpen()
    self.listWidget = None
    if (len(self.nameList)>1):
      self.listWidget = QtGui.QListWidget(self)
      for name in self.nameList:
        self.listWidget.addItem(name)
      #self.listWidget.setSelectionMode(
      #    QtGui.QListWidget.ExtendedSelection)

    # Set up the central tabbed widget.
    tabBar = QtGui.QTabWidget(self)
    inputScreen = InputWidget(self)
    tabBar.addTab(inputScreen,"Input")
    execScreen = QtGui.QFrame(self)
    tabBar.addTab(execScreen,"Execute")
    outputScreen = OutputWidget(self)
    tabBar.addTab(outputScreen,"Output")

    # Connect ListView to InputWidget and OutputWidget.
    if self.listWidget:
      QtCore.QObject.connect(self.listWidget.selectionModel(),
          QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),
          inputScreen, QtCore.SLOT("respondToSelection()"))
      QtCore.QObject.connect(self.listWidget.selectionModel(),
          QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),
          outputScreen, QtCore.SLOT("respondToSelection()"))

    hbox = QtGui.QHBoxLayout()
    if self.listWidget != None:
      splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
      splitter.addWidget(self.listWidget)
      splitter.addWidget(tabBar)
      splitter.setSizes([150,850])
      hbox.addWidget(splitter)
      self.listWidget.item(0).setSelected(True)
    else:
      hbox.addWidget(tabBar)
      inputScreen.readXML()
      outputScreen.readH5()
    self.setLayout(hbox)

  def autoOpen(self):
    self.nameList = glob.glob("pimc.xml")
    if len(self.nameList)==1:
       self.nameList[0] = "./"
    else: 
      for file in glob.glob("*/pimc.xml"):
        self.nameList.append(file[:-8])
      for file in glob.glob("*/*/pimc.xml"):
        self.nameList.append(file[:-8])
      for file in glob.glob("*/*/*/pimc.xml"):
        self.nameList.append(file[:-8])
      for file in glob.glob("*/*/*/*/pimc.xml"):
        self.nameList.append(file[:-8])
      for file in glob.glob("*/*/*/*/*/pimc.xml"):
        self.nameList.append(file[:-8])
      for file in glob.glob("*/*/*/*/*/*/pimc.xml"):
        self.nameList.append(file[:-8])
      self.nameList.sort()
