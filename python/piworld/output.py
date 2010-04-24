import sys,os,math,numpy
import pitools
from pitools import Unit
from PyQt4 import QtGui,QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from estimatorview import *
from conductivityview import *
from scalarview import *
from densityview import *
from permutationview import *

class SimulationData:
  def __init__(self):
    self.temperature = 0.
    self.file = None
    self.estimators = []

class OutputWidget(QtGui.QWidget):
  def __init__(self, parent):
    QtGui.QWidget.__init__(self,parent)
    self.project = parent
    self.project.oldSelectedEstimatorNames = []
    self.data = [None]*len(self.project.nameList)

    self.setContentsMargins(2,2,2,2)
    vbox = QtGui.QVBoxLayout()
    vbox.setContentsMargins(0,0,0,0)
    vbox.setSpacing(3)
    self.simInfoWidget = SimInfoWidget(self)
    vbox.addWidget(self.simInfoWidget,1)
    vbox.addStretch(1)
    self.estimatorWidget = EstimatorWidget(self.project,self)
    vbox.addWidget(self.estimatorWidget,10)
    self.setLayout(vbox)

  def readH5(self,index=0):
    if self.data[index] == None:
      try:
        dat = SimulationData()
        dat.file = pitools.openFile(self.project.nameList[index]+"pimc.h5")
        dat.temperature = dat.file.getTemperature(Unit.K) 
        dat.nslice = dat.file.getNSlice() 
        dat.superCell = dat.file.getSuperCell(Unit.nm) 
        dat.estimators = dat.file.getEstimatorList()
        #dat.file.close()
        self.data[index] = dat
      except IOError:
        pass
    self.myupdate(self.data[index])

  @QtCore.pyqtSlot()
  def respondToSelection(self):
    selection = map(QtCore.QModelIndex.row,
                    self.project.listWidget.selectedIndexes())
    selection.sort()
    if len(selection) > 0:
      self.readH5(selection[0])
    else:
      pass

  def myupdate(self,data):
    self.simInfoWidget.myupdate(data)
    self.estimatorWidget.myupdate(data)

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
  def __init__(self, project, parent=None):
    QtGui.QWidget.__init__(self,parent)
    self.project = project
    self.layout = QtGui.QHBoxLayout()
    self.layout.setSpacing(1)
    self.setContentsMargins(1,1,1,1)
    splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
    self.listWidget = QtGui.QListWidget(parent) 
    self.listWidget.setContentsMargins(1,1,1,1)
    splitter.addWidget(self.listWidget)
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
    oldSelectedNames = [item.text() for item in self.listWidget.selectedItems()]
    if len(oldSelectedNames)==0:
      oldSelectedNames = self.project.oldSelectedEstimatorNames
    else:
      self.project.oldSelectedEstimatorNames = oldSelectedNames
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
    for name in oldSelectedNames:
      for item in self.listWidget.findItems(name, QtCore.Qt.MatchExactly):
        item.setSelected(True)

  @QtCore.pyqtSlot()
  def respondToSelection(self):
    selection = map(QtCore.QModelIndex.row,
                    self.listWidget.selectedIndexes())
    selection.sort()
    if len(selection) > 0:
      name = self.listWidget.item(selection[0]).text()
      est = self.getEstimatorViewer(name)
      if not getattr(est,"vindex",False):
        if est.view == None:
          est.vindex = 0
        else:
          est.vindex = self.estimatorView.count()
          self.estimatorView.addWidget(est.view)
      self.estimatorView.setCurrentIndex(est.vindex)
    else:
      self.estimatorView.setCurrentIndex(0)

  def getEstimatorViewer(self, name):
    est = (e for e in self.data.estimators if e.name == name).next() 
    if getattr(est,"view","missing") == "missing":
      est.view = None
      if est.type == 258:
        est.view = ConductivityView(est,self.data)
      elif est.type >= 64 and est.type <128:
        est.view = ScalarView(est,self.data)
      elif est.type == 129:
        est.view = DensityView(est,self.data)
      elif est.type == 1:
        est.view = PermutationView(est,self.data)
      else:
        print "Mising viewer type %i for %s." % (est.type,est.name)
    return est
