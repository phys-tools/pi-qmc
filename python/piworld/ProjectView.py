import sys,os,glob
from PyQt4 import QtGui,QtCore
from input import InputWidget
from output import OutputWidget

class ProjectView(QtGui.QWidget):
  def __init__(self, projectModel, parent=None):
    QtGui.QWidget.__init__(self,parent)
    self.model = projectModel
    self.listWidget = None
    if (self.model.hasMultipleSimulations):
      self.listWidget = QtGui.QListWidget(self)
      for name in self.model.nameList:
        self.listWidget.addItem(name)
      #self.listWidget.setSelectionMode(
      #    QtGui.QListWidget.ExtendedSelection)

    # Set up the central tabbed widget.
    tabBar = QtGui.QTabWidget(self)
    inputScreen = InputWidget(self,self.model)
    tabBar.addTab(inputScreen,"Input")
    execScreen = QtGui.QFrame(self)
    tabBar.addTab(execScreen,"Execute")
    outputScreen = OutputWidget(self,self.model)
    tabBar.addTab(outputScreen,"Output")

    # Connect ListView to ProjectModel.
    if self.listWidget:
      QtCore.QObject.connect(self.listWidget.selectionModel(),
          QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),
          self.model, QtCore.SLOT("respondToSelection(QItemSelection)"))

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
      inputScreen.respondToSelection()
      outputScreen.respondToSelection()
    self.setLayout(hbox)
