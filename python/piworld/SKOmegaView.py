#-*- coding: utf-8 -*-
import math,scipy,numpy as np,tables,pylab
import math,scipy,numpy as np,tables,matplotlib.pyplot as plt
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *


class SKOmegaView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    #Read the data and do some setup.
    self.polk = estimatorNode.node.read()
    self.boxl = data.superCell
    #Divide correlation by total volume to get polarizability (definition).
    self.polk /= self.boxl.prod()
    # Make grid of k values for plotting and analysis.
    # Uses some python magic, especialy map and ogrid.
    nx = self.polk.shape[2:-1]
    def myslice(n): return slice(-n/2,n/2)
    kval = (2*math.pi/self.boxl)*np.array(np.ogrid[map(myslice,nx)])
    self.kgrid = np.sqrt(np.fft.fftshift((kval**2).sum(0)))

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
      #Pick off dimensions of the grid.
      nx = self.data.polk.shape[2:-1]
      npart = self.data.polk.shape[0]
      self.axes.axhline(0)
      ifreq = 0
      kaxis = self.data.kgrid.reshape(-1)
      for ipart in xrange(npart):
        for jpart in xrange(npart):
          intra = np.real(self.data.polk[ipart,jpart,...,ifreq])
          self.axes.plot(kaxis[1:], intra.reshape(-1)[1:], ls='None', 
                         marker=".",ms=1.5)
      self.axes.set_xlabel(r"$q$")
      self.axes.set_ylabel(r"$\chi_{nn}$")
      #self.axes.axis(xmin=1e-9)
