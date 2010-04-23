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

    pylab.rc('font', family='serif', size=9)

    self.setContentsMargins(1,1,1,1)

    self.tabBar = QtGui.QTabWidget(self)
    self.tabBar.setTabPosition(QtGui.QTabWidget.South)

    self.tracePlot =  self.TracePlotWidget(self)
    self.tabBar.addTab(self.tracePlot,"Trace")
    self.autocorPlot =  self.PlotWidget(self)
    self.tabBar.addTab(self.autocorPlot,"Autocorrelation")
    self.blockingPlot =  self.PlotWidget(self)
    self.tabBar.addTab(self.blockingPlot,"Blocking")

    hbox = QtGui.QHBoxLayout()
    hbox.setSpacing(0)
    hbox.addWidget(self.tabBar,15)

    self.info = self.InfoWidget(self)
    hbox.addWidget(self.info,1)

    self.setLayout(hbox)
    

  class TracePlotWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      values = self.data.scalar.data
      times = numpy.arange(len(values))
      # Compute average and error, and plot trace.
      av,err = self.data.scalar.getAverage()
      ax = self.axes = self.figure.add_axes([0.15,0.15,0.72,0.82])
      ax.axhspan(av-err,av+err,color='y',alpha=0.2)
      ax.axhline(y=av,color='k',lw=0.2)
      ax.plot(times,values)
      ax.axis(xmax=len(values)-1)
      ax.set_xlabel(r"Monte Carlo timestep")
      ax.set_ylabel(self.data.scalar.name)
      # Compute histogram.
      hax = self.histAxes \
          = self.figure.add_axes([0.87,0.15,0.12,0.82],xticks=[],yticks=[])
      hax.hist(values,bins=32,orientation='horizontal',
               range=ax.get_ybound())
       

  class PlotWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      pass

  class InfoWidget(QtGui.QWidget): 
    def __init__(self, data, parent=None):
      self.data = data
      QtGui.QWidget.__init__(self,parent)
 
      mean,error = self.data.scalar.getAverage()

      views = []
      labels = []

      self.meanView = QtGui.QLabel("%g"%mean)
      labels.append(u"mean"); views.append(self.meanView)
      self.errorView = QtGui.QLabel("%g"%error)
      labels.append(u"error of mean"); views.append(self.errorView)
      self.sigmaView = QtGui.QLabel("???")
      labels.append(u"σ"); views.append(self.sigmaView)
      self.autocorView = QtGui.QLabel("???")
      labels.append(u"T<sub>autocorr.</sub>"); views.append(self.autocorView)
      self.startView = QtGui.QLabel("0")
      labels.append(u"start cutoff"); views.append(self.startView)
      self.endView = QtGui.QLabel("?")
      labels.append(u"end cutoff"); views.append(self.endView)

      vbox = QtGui.QVBoxLayout()
      vbox.setSpacing(1)

      for view,label in zip(views,labels):
        vbox.addWidget(QtGui.QLabel(u"<font size='-1'>%s</font>"%label))
        vbox.addWidget(view)
        vbox.addStretch(1)

      vbox.addStretch(1)
      self.setLayout(vbox)
