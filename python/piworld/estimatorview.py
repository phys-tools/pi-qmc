import sys,os,math,numpy
from PyQt4 import QtGui,QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class EstimatorView(QtGui.QWidget):
  def __init__(self, parent=None):
    QtGui.QWidget.__init__(self,parent)


class MyMplCanvas(FigureCanvasQTAgg):
  def __init__(self, parent=None, width=5, height=4, dpi=100):
    self.figure = Figure(figsize=(width,height), dpi=dpi)
    self.figure.set_facecolor('w')
    self.computeInitialFigure()  
    FigureCanvasQTAgg.__init__(self,self.figure)
    self.setParent(parent)
  def computeInitialFigure(self):
    pass

class MyGraph(MyMplCanvas):
  def computeInitialFigure(self):
    t = numpy.arange(0.0,3.0,0.01)
    s = numpy.sin(2*math.pi*t)
    self.axes = self.figure.set_axes([0.2,0.2,0.79,0.79])
    self.axes.plot(t,s)
