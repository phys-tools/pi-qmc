import sys,os,glob
from PyQt4 import QtGui,QtCore
from SimulationModel import *
from lxml import etree

class ProjectModel(QtCore.QObject):
  def __init__(self,gui):
    QtCore.QObject.__init__(self,None)    
    self.autoOpen()
    self.prefs = {}   
    self.gui = gui
  
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
      # Check for remote access info in pimc.xml.
      xmldoc = etree.fromstring(self.getInputXMLText())
      hostNode = xmldoc.xpath("Host")
      hostMachine = None
      hostDirectory = None
      if len(hostNode)>0:
        hostMachine=hostNode[0].attrib.get('machine')
        hostDirectory=hostNode[0].attrib.get('directory')
      # Load the data.
      self.data[index] = dataFromH5File(self,self.nameList[index]+"pimc.h5")
      self.data[index].xmldoc = xmldoc
      if hostMachine: self.data[index].hostMachine = hostMachine
      if hostDirectory: self.data[index].hostDirectory = hostDirectory
      QtCore.QObject.connect(self.data[index],
        QtCore.SIGNAL('statusMessage(QString,int)'),
        self.gui.statusbar, QtCore.SLOT('message(QString,int)'))
    return self.data[index]

  @QtCore.pyqtSlot()
  def respondToDataChange(self):
    self.selectionChanged.emit()

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
    if len(self.nameList)>0: self.selection=[0]
    self.xmlstring = [None]*len(self.nameList)
    self.data = [None]*len(self.nameList)
