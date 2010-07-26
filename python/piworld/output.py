import sys,os,math,numpy
import pitools
from pitools import Unit
from PyQt4 import QtGui,QtCore
from EstimatorView import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class OutputWidget(QtGui.QWidget):
  def __init__(self, parent, projectModel):
    QtGui.QWidget.__init__(self,parent)
    self.projectModel = projectModel
    self.projectModel.prefs['LastSelectedEstimatorNames'] = []
    self.currentIndex = 0

    self.setContentsMargins(2,2,2,2)
    vbox = QtGui.QVBoxLayout()
    vbox.setContentsMargins(0,0,0,0)
    vbox.setSpacing(3)
    self.simInfoWidget = SimInfoWidget(self)
    vbox.addWidget(self.simInfoWidget,1)
    vbox.addStretch(1)
    self.estimatorWidget = EstimatorWidget(self.projectModel,self)
    vbox.addWidget(self.estimatorWidget,10)
    self.setLayout(vbox)
    # Watch model selection.
    QtCore.QObject.connect(self.projectModel,
        QtCore.SIGNAL("selectionChanged()"),
        self, QtCore.SLOT("respondToSelection()"))

  @QtCore.pyqtSlot()
  def respondToSelection(self):
    print "Output responding to selection"    
    self.data = self.projectModel.getSimulationData()
    self.simInfoWidget.myupdate(self.data)
    self.estimatorWidget.myupdate(self.data)


class SimInfoWidget(QtGui.QGroupBox):
  def __init__(self, parent=None):
    QtGui.QWidget.__init__(self,parent)
    self.setTitle("Simulation Info")
    hbox = QtGui.QHBoxLayout()
    hbox.setSpacing(3)
    self.labelT = QtGui.QLabel("T=? K (? slices)",self)
    hbox.addWidget(self.labelT)
    hbox.addStretch(1)
    self.labelSuperCell = QtGui.QLabel("Box Size: (?)",self)
    hbox.addWidget(self.labelSuperCell)
    hbox.addStretch(10)
    self.setLayout(hbox)

  def myupdate(self,data):
    self.data = data
    if data==None:
      self.labelT.setText("T=? (? slices)")
      self.labelSuperCell.setText("Box Size: (?)")
    else:
      self.labelT.setText("T=%g K (%i slices)"
          % (data.temperature, data.nslice))
      sctext = "Box Size: ( %g" % data.superCell[0]
      for x in data.superCell[1:]:
        sctext += (", %g" % x)
      self.labelSuperCell.setText(sctext+" ) nm")

class EstimatorWidget(QtGui.QWidget):
  def __init__(self, projectModel, parent=None):
    QtGui.QWidget.__init__(self,parent)
    self.project = projectModel
    self.layout = QtGui.QHBoxLayout()
    self.layout.setSpacing(1)
    self.setContentsMargins(1,1,1,1)
    splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
    self.listWidget = QtGui.QListWidget(self) 
    self.listWidget.setContentsMargins(1,1,1,1)
    panel = QtGui.QWidget()
    vbox = QtGui.QVBoxLayout()
    panel.setLayout(vbox)
    vbox.addWidget(self.listWidget)
    splitter.addWidget(panel)
    self.estimatorView = QtGui.QStackedWidget(self)
    splitter.addWidget(self.estimatorView)
    self.layout.addWidget(splitter)
    self.plot = EstimatorView(self)
    self.estimatorView.addWidget(self.plot)
    self.setLayout(self.layout)
    splitter.setSizes([200,800])
    QtCore.QObject.connect(self.listWidget.selectionModel(),
        QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),
        self, QtCore.SLOT("respondToSelection()"))

  def myupdate(self,data):
    # First save name of previously selected estimator.
    oldNames  = map(QtGui.QListWidgetItem.text, self.listWidget.selectedItems())
    if len(oldNames)==0:
      oldNames = self.project.prefs['LastSelectedEstimatorNames']
    else:
      self.project.prefs['LastSelectedEstimatorNames'] = oldNames
    # Now remove old list items are replace with current data.
    self.listWidget.setSelectionMode(QtGui.QListWidget.NoSelection)
    self.data = data
    if data==None:
      if self.listWidget.count() > 0:
        for i in range(self.listWidget.count()-1,-1,-1):
          self.listWidget.takeItem(i)
    else:
      if self.listWidget.count() > 0:
        for i in range(self.listWidget.count()-1,-1,-1):
          self.listWidget.takeItem(i)
      for node in data.estimators:
        if node.type >= 0:
          self.listWidget.addItem(node.name)
    # If name matches previously selected estimator, select it.
    self.listWidget.setSelectionMode(QtGui.QListWidget.SingleSelection)
    for name in oldNames:
      for item in self.listWidget.findItems(name, QtCore.Qt.MatchExactly):
        item.setSelected(True)

  @QtCore.pyqtSlot()
  def respondToSelection(self):
    selection = map(QtCore.QModelIndex.row,
                    self.listWidget.selectedIndexes())
    selection.sort()
    if len(selection) > 0:
      name = self.listWidget.item(selection[0]).text()
      est = self.data.getEstimatorViewer(name)
      if not getattr(est,"vindex",False):
        if est.view == None:
          est.vindex = 0
        else:
          est.vindex = self.estimatorView.count()
          self.estimatorView.addWidget(est.view)
      self.estimatorView.setCurrentIndex(est.vindex)
    else:
      self.estimatorView.setCurrentIndex(0)
