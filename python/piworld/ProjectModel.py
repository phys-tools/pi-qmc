import sys,os,glob
from PyQt4 import QtGui,QtCore
from SimulationModel import *

class ProjectModel(QtCore.QObject):
  def __init__(self):
    QtCore.QObject.__init__(self,None)    
    self.autoOpen()
    self.prefs = {}   
  
  currentIndex = 0
  selection = []

  # Signals
  selectionChanged = QtCore.pyqtSignal()

  # Slots
  @QtCore.pyqtSlot(QtGui.QItemSelection)
  def respondToSelection(self, selection):
    self.selection = map(QtCore.QModelIndex.row, selection.indexes())
    self.selection.sort()
    self.selectionChanged.emit()

  def getInputXMLText(self):
    index = 0
    if self.hasMultipleSimulations:
      if len(self.selection) == 0: return ""
      index = self.selection[0]
    if self.xmlstring[index] == None:
      file = open(self.nameList[index]+"pimc.xml","r")    
      self.xmlstring[index] = file.read()
      file.close()
    return self.xmlstring[index]

  def getSimulationData(self):
    index = 0
    if self.hasMultipleSimulations:
      if len(self.selection) == 0: return None
      index = self.selection[0]
    if self.data[index] == None:
      self.data[index] = dataFromH5File(self,self.nameList[index]+"pimc.h5")
    return self.data[index]

  def autoOpen(self):
    self.nameList = glob.glob("pimc.xml")
    if len(self.nameList)==1:
      self.hasMultipleSimulations = False
      self.nameList[0] = "./"
    else: 
      self.hasMultipleSimulations = True
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
    self.xmlstring = [None]*len(self.nameList)
    self.data = [None]*len(self.nameList)
