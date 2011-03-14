#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from PyQt4 import QtCore, QtGui
from piworldUI import Ui_piworld
from ProjectModel import ProjectModel
from MyHighlighter import MyHighlighter

class StartQT4(QtGui.QMainWindow):
  def __init__(self, parent=None):
    QtGui.QWidget.__init__(self, parent)
    self.ui = Ui_piworld()
    self.ui.setupUi(self)
    self.statusBar().showMessage("Initializing...",3000)

    # Cosmetic adjustments to gui.
    self.ui.mainTabWidget.setCurrentIndex(0)
    #self.ui.mainTabWidget.setTabEnabled(1,False)
    self.ui.estimatorSplitter.setSizes([150,850])
    self.ui.outputSplitter.setSizes([200,800])
    self.ui.simulationList.hide()
    self.ui.menuView.addAction(self.ui.toolBar.toggleViewAction())
    self.ui.toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)
    
    # Set up the project model.
    self.project = ProjectModel(self.ui)
    self.project.prefs['LastSelectedEstimatorNames'] = []

    # Set up the simulation list and connect to ProjectModel
    if (self.project.hasMultipleSimulations):
      for name in self.project.nameList:
        self.ui.simulationList.addItem(name)
      self.ui.simulationList.show()
      self.ui.simulationList.setEnabled(True)
      self.ui.simulationList.setCurrentRow(0)
      QtCore.QObject.connect(self.ui.simulationList.selectionModel(),
        QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),
        self.project, QtCore.SLOT("respondToSelection(QItemSelection)"))

    # Watch project model selection.
    QtCore.QObject.connect(self.project, QtCore.SIGNAL("selectionChanged()"),
        self, QtCore.SLOT("respondToSimulationSelection()"))

    # Watch estimator selection.
    QtCore.QObject.connect(self.ui.estimatorList,
      QtCore.SIGNAL("itemSelectionChanged()"),
      self, QtCore.SLOT("respondToEstimatorSelection()"))

    # Watch reload.
    QtCore.QObject.connect(self.ui.actionRefresh,
      QtCore.SIGNAL("triggered(bool)"),
      self, QtCore.SLOT("respondToRefresh()"))

    # Set up the input text.
    MyHighlighter(self.ui.inputXMLText)
    self.ui.inputXMLText.setPlainText(self.project.getInputXMLText())

    self.project.selectionChanged.emit()

  @QtCore.pyqtSlot()
  def respondToSimulationSelection(self):
    #
    # Update input text.
    #
    self.ui.inputXMLText.setPlainText(self.project.getInputXMLText())
    #
    # Update the current simulation data settings.
    #
    data = self.project.getSimulationData()
    self.ui.actionRefresh.setEnabled(True)
    #
    # Update simulation info.
    #
    if data==None:
      self.ui.temperatureField.setText("?")
      self.ui.superCellField.setText("(?)")
      self.ui.nslicesField.setText("?")
    else:
      self.ui.temperatureField.setText("%g K"%data.temperature)
      sctext = "( %g" % data.superCell[0]
      for x in data.superCell[1:]:
        sctext += (", %g" % x)
      self.ui.superCellField.setText(sctext+" ) nm")
      self.ui.nslicesField.setText("%i"%data.nslice)
    #
    # Update the estimator browser.
    #
    # First save name of previously selected estimators.
    oldNames  = map(QtGui.QListWidgetItem.text,
                    self.ui.estimatorList.selectedItems())
    if len(oldNames)==0:
      oldNames = self.project.prefs['LastSelectedEstimatorNames']
    else:
      self.project.prefs['LastSelectedEstimatorNames'] = oldNames
    # Now remove old list items and replace with current data.
    self.ui.estimatorList.setSelectionMode(QtGui.QListWidget.NoSelection) 
    if self.ui.estimatorList.count() > 0:
      for i in range(self.ui.estimatorList.count()-1,-1,-1):
        self.ui.estimatorList.takeItem(i)
    if data:
      for node in data.estimators:
        if node.typeString[-4:] != "-err":
          self.ui.estimatorList.addItem(node.name)
    # Remove estimator views from the estimatorViewStack.
    while self.ui.estimatorViewStack.count()>2:
      self.ui.estimatorViewStack.removeWidget(
        self.ui.estimatorViewStack.widget(
          self.ui.estimatorViewStack.count()-1))
    # If name matches previously selected estimator, select it.
    self.ui.estimatorList.setSelectionMode(QtGui.QListWidget.SingleSelection)
    for name in oldNames:
      for item in self.ui.estimatorList.findItems(name, QtCore.Qt.MatchExactly):
        item.setSelected(True)

  @QtCore.pyqtSlot()
  def respondToEstimatorSelection(self):
    selection = self.ui.estimatorList.selectedItems()
    if len(selection)==0:
      self.ui.estimatorViewStack.setCurrentWidget(self.ui.emptyEstimatorView)
      return
    name = selection[0].text()
    est = self.project.getSimulationData().getEstimatorViewer(name)
    if est.view == None:
      self.ui.estimatorViewStack.setCurrentWidget(self.ui.missingEstimatorView)
    else:
      index = self.ui.estimatorViewStack.indexOf(est.view)
      if index==-1:
        self.ui.estimatorViewStack.addWidget(est.view)
        index = self.ui.estimatorViewStack.indexOf(est.view)
      self.ui.estimatorViewStack.setCurrentIndex(index)

  @QtCore.pyqtSlot()
  def respondToRefresh(self):
    self.statusBar().showMessage("Reloading data...",1000)
    self.project.getSimulationData().reloadData()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = StartQT4()
    myapp.show()
    sys.exit(app.exec_())
