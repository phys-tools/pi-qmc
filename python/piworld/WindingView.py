#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
import math,scipy,numpy,tables,matplotlib.pyplot as plt
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *


class WindingView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    self.histogram = estimatorNode.node.read()

    self.plot =  self.PlotWidget(self)
    vbox = QtGui.QVBoxLayout()
    vbox.setSpacing(0)
    vbox.addWidget(self.plot,10)

    self.setLayout(vbox)
    
  class PlotWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.14,0.15,0.69,0.83])
      plt.rc('font', family='serif', size=9)
      self.cbaxes = self.figure.add_axes([0.84,0.15,0.05,0.83])
      nmax = int(self.data.histogram.shape[0]/2)
      im =  self.axes.imshow(self.data.histogram,
              extent=[-nmax,nmax,-nmax,nmax],
              origin="lower", interpolation="nearest",
	      vmin=0.0)
      plt.colorbar(im,cax=self.cbaxes)
