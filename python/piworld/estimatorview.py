import sys,os,math,numpy
from PyQt4 import QtGui,QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class EstimatorView(QtGui.QWidget):
  def __init__(self, parent=None):
    QtGui.QWidget.__init__(self,parent)


class MyMplCanvas(FigureCanvasQTAgg):
  def __init__(self, parent=None, width=5, height=4, dpi=100):
    figure = Figure(figsize=(width,height), dpi=dpi)
    self.axes = figure.add_subplot(111)
    self.axes.hold(True)
    self.computeInitialFigure()  
    FigureCanvasQTAgg.__init__(self,figure)
    self.setParent(parent)
  def computeInitialFigure(self):
    pass

class MyGraph(MyMplCanvas):
  def computeInitialFigure(self):
    t = numpy.arange(0.0,3.0,0.01)
    s = numpy.sin(2*math.pi*t)
    self.axes.plot(t,s)

