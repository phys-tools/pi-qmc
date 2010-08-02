import sys,os,glob
from PyQt4 import QtGui,QtCore
import pitools
from pitools import Unit
from EstimatorView import *
from ConductivityView import *
from ScalarView import *
from DensityView import *
from PermutationView import *
from SimulationModel import *

class SimulationModel(QtCore.QObject):
  def __init__(self, project):
    QtCore.QObject.__init__(self,None)
    self.project = project
    self.temperature = 0.
    self.nslice = 0
    self.estimators = []
    self.superCell = None
    self.file = None
    QtCore.QObject.connect(self, QtCore.SIGNAL("dataChanged()"),
      self.project, QtCore.SLOT("respondToDataChange()"))
    

  # Signals
  dataChanged = QtCore.pyqtSignal()

  #Slots
  @QtCore.pyqtSlot()
  def reloadData(self):
    self.file.close()
    try:
      self.file = pitools.openFile(self.filename)
      self.temperature = self.file.getTemperature(Unit.K) 
      self.nslice = self.file.getNSlice() 
      self.superCell = self.file.getSuperCell(Unit.nm) 
      self.estimators = self.file.getEstimatorList()
    except IOError:
      return None
    self.dataChanged.emit()

  def getEstimatorViewer(self, name):
    est = (e for e in self.estimators if e.name == name).next() 
    if getattr(est,"view","missing") == "missing":
      est.view = None
      if est.type == 258:
        est.view = ConductivityView(est,self)
      elif est.type >= 64 and est.type <128:
        est.view = ScalarView(est,self)
      elif est.type == 129:
        est.view = DensityView(est,self)
      elif est.type == 1:
        est.view = PermutationView(est,self)
      else:
        print "Mising viewer type %i for %s." % (est.type,est.name)
    return est

def dataFromH5File(project, filename):
  data = SimulationModel(project)
  data.filename = filename
  try:
    data.file = pitools.openFile(data.filename)
    data.temperature = data.file.getTemperature(Unit.K) 
    data.nslice = data.file.getNSlice() 
    data.superCell = data.file.getSuperCell(Unit.nm) 
    data.estimators = data.file.getEstimatorList()
  except IOError:
    return None
  return data

