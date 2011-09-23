import sys,os,glob
from PyQt4 import QtGui,QtCore
import pitools
from pitools import Unit
from EstimatorView import *
from ConductivityView import *
from ScalarView import *
from DensityView import *
from PairCFView import *
from FreeEnergyView import *
from WindingView import *
from SKOmegaView import *
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
  # Signals
  statusMessage = QtCore.pyqtSignal(QtCore.QString,int)

  #Slots
  @QtCore.pyqtSlot()
  def reloadData(self):
    self.file.close()
    # Retreive remote data if needed.
    if hasattr(self,"hostMachine") and hasattr(self,"hostDirectory"):
      self.statusMessage.emit("downloaded new data from %s:%s"
        % (self.hostMachine, self.hostDirectory),5000)
      os.system('scp "%s:%s/pimc.h5" %s'
                 % (self.hostMachine, self.hostDirectory,
                    self.project.nameList[self.index]))
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
      if est.typeString == "dynamic-array/conductivity":
        est.view = ConductivityView(est,self)
      if est.typeString == "dynamic-array/structure-factor":
        est.view = SKOmegaView(est,self)
      elif est.typeString[:6] == "scalar":
        est.view = ScalarView(est,self)
      elif est.typeString == "array/density":
        est.view = DensityView(est,self)
      elif est.typeString == "histogram/permutation":
        est.view = PermutationView(est,self)
      elif est.typeString == "histogram/free-energy":
        est.view = FreeEnergyView(est,self)
      elif est.typeString == "histogram/winding":
        est.view = WindingView(est,self)
      elif est.typeString == "array/pair-correlation":
        est.view = PairCFView(est,self)
      else:
        print "Missing viewer type %s for %s." % (est.typeString,est.name)
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
