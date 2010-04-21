#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from estimatorview import *


class ScalarView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    self.scalar = data.file.getScalar(estimatorNode.name) 
    #print self.scalar.getAverage()
    #print self.scalar.unit

    pylab.rc('font', family='serif', size=9)

    self.tabBar = QtGui.QTabWidget(self)
    self.tabBar.setTabPosition(QtGui.QTabWidget.South)

    self.tracePlot =  self.PlotWidget(self)
    self.tabBar.addTab(self.tracePlot,"Trace")
    self.autocorPlot =  self.PlotWidget(self)
    self.tabBar.addTab(self.autocorPlot,"Autocorrelation")
    self.blockingPlot =  self.PlotWidget(self)
    self.tabBar.addTab(self.blockingPlot,"Blocking")

    hbox = QtGui.QHBoxLayout()
    hbox.setSpacing(0)
    hbox.addWidget(self.tabBar,10)

    self.info = self.InfoWidget(self)
    hbox.addWidget(self.info,1)

    self.setLayout(hbox)
    

  class PlotWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      pass
      #maxfreq = math.ceil((self.data.omegan/self.data.omega1).max())
      #maxg = math.ceil((self.data.sigma0iw).max())
      #omegas = numpy.arange(0,maxfreq+0.4,0.1)*self.data.omega1
      #self.axes.errorbar(self.data.omegan[:]/self.data.omega1, 
      #                   self.data.sigma0iw[self.data.ndx/2,:],
      #                   self.data.sigma0iw_err[self.data.ndx/2,:], fmt=".")
      #self.axes.plot(omegas/self.data.omega1,
      #               self.data.peval(omegas,self.data.pfit_a))
      #self.axes.axis([0,maxfreq+0.4,0,maxg+0.4])
      #self.axes.set_xlabel(r"$\omega_n/\omega_1$")
      #self.axes.set_ylabel(r"$G$ ($e^2/h$)")

  class InfoWidget(QtGui.QWidget): 
    def __init__(self, data, parent=None):
      self.data = data
      QtGui.QWidget.__init__(self,parent)
 
      #gfit = self.data.peval(0.,self.data.pfit_a)

      vbox = QtGui.QVBoxLayout()
      vbox.setSpacing(0)
      #self.labelMean = QtGui.QLabel(u"mean")
      vbox.addWidget(QtGui.QLabel(u"mean"))
      mean,error = self.data.scalar.getAverage()
      self.meanView = QtGui.QLabel("%g"%mean)
      vbox.addWidget(self.meanView)
      vbox.addWidget(QtGui.QLabel(u"error of mean"))
      self.errorView = QtGui.QLabel("%g"%error)
      vbox.addWidget(self.errorView)
      vbox.addWidget(QtGui.QLabel(u"σ"))
      vbox.addWidget(QtGui.QLabel(u"autocor. time"))
      vbox.addWidget(QtGui.QLabel(u"start cutoff"))
      vbox.addWidget(QtGui.QLabel(u"end cutoff"))
      vbox.addStretch(1)
      self.setLayout(vbox)
