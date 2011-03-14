#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
import math,scipy,numpy,tables,matplotlib.pyplot as plt
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *


class FreeEnergyView(EstimatorView):

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
      self.axes = self.figure.add_axes([0.14,0.15,0.84,0.83])
      plt.rc('font', family='serif', size=9)
      nbin = len(self.data.histogram)
      x = numpy.arange(nbin)
      prob = self.data.histogram/self.data.histogram.sum()
      firstProb = self.data.histogram[0] / \
                 (self.data.histogram*numpy.arange(1,nbin+1)).sum()
      self.axes.bar(x-0.4, prob)
      self.axes.set_xlabel(r"model")
      self.axes.set_ylabel(r"$\exp(-\Delta F/k_B T)$")
      self.axes.axis(xmin=-0.495, xmax=nbin-0.505,
                     ymin=0., ymax=numpy.max(prob)*1.1)
