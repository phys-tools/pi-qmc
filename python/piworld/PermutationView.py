#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *


class PermutationView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    self.permutation = estimatorNode.node.read()

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
      self.axes = self.figure.add_axes([0.18,0.15,0.80,0.83])
      nbin = len(self.data.permutation)
      x = numpy.arange(nbin)+1
      prob = self.data.permutation/self.data.permutation.sum()
      self.axes.bar(x-0.4, prob)
      self.axes.set_xlabel(r"cycle length")
      self.axes.set_ylabel(r"probability")
      self.axes.axis(xmin=0, xmax=nbin+0.495, ymax=numpy.max(prob)*1.1)
