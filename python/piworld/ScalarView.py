#-*- coding: utf-8 -*-
from PyQt4 import QtCore, QtGui
import math,scipy,numpy,tables,matplotlib.pyplot as plt
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *
from ScalarViewUI import Ui_ScalarView



class ScalarView(QtGui.QWidget):
  def __init__(self, estimatorNode, data, parent=None):
    QtGui.QWidget.__init__(self,parent)
    self.ui = Ui_ScalarView()
    self.ui.setupUi(self)
    # Read scalar data.
    self.scalar = data.file.getScalar(estimatorNode.name) 
    self.model = self.ScalarAnalysisModel(self.scalar)
    # Set up plots in tabs.
    self.tracePlot =  self.TracePlotWidget(self.model)
    self.ui.tabWidget.addTab(self.tracePlot,"Trace")
    self.autocorrPlot =  self.AutocorrPlotWidget(self.model)
    self.ui.tabWidget.addTab(self.autocorrPlot,"Autocorrelation")
    self.blockingPlot =  self.BlockingPlotWidget(self.model)
    self.ui.tabWidget.addTab(self.blockingPlot,"Blocking")
    # Make connections and update data fields.
    QtCore.QObject.connect(self.ui.startText,
      QtCore.SIGNAL("editingFinished()"),
      self, QtCore.SLOT("respondToStartChange()"))
    QtCore.QObject.connect(self.ui.endText,
      QtCore.SIGNAL("editingFinished()"),
      self, QtCore.SLOT("respondToEndChange()"))
    self.updateFields()
    self.ui.splitter.setSizes([920,80])

  def updateFields(self):
    self.ui.meanText.setText("%g"%self.model.av)
    self.ui.errorText.setText("%g"%self.model.err)
    self.ui.sigmaText.setText("%g"%self.model.sigma)
    self.ui.autocorrText.setText("%g"%self.model.actime)
    self.ui.startText.setText("%g"%self.model.start)
    self.ui.endText.setText("%g"%self.model.end)

  @QtCore.pyqtSlot()
  def respondToStartChange(self):
    oldStart = self.model.start
    start = int(self.ui.startText.text())
    if start >= self.model.end: start = self.model.end-1
    if start < 0: start = 0
    self.ui.startText.setText("%i"%start)
    if start != oldStart:
      self.model.start = start
      self.model.calculate()
      self.respondToModelChange()

  @QtCore.pyqtSlot()
  def respondToEndChange(self):
    oldEnd = self.model.end
    end = int(self.ui.endText.text())
    if end <= self.model.start: end = self.model.start+1
    if end >= self.model.ndata: end = self.model.ndata-1
    self.ui.endText.setText("%i"%end)
    if end != oldEnd:
      self.model.end = end
      self.model.calculate()
      self.respondToModelChange()

  @QtCore.pyqtSlot()
  def respondToModelChange(self):
    self.updateFields()
    self.tracePlot.plotFigure(); self.tracePlot.draw()
    self.autocorrPlot.plotFigure(); self.autocorrPlot.draw()
    self.blockingPlot.plotFigure(); self.blockingPlot.draw()

  class ScalarAnalysisModel:
    def __init__(self, scalar):
      self.scalar = scalar
      self.values = self.scalar.data
      self.ndata = len(self.values)
      self.start, self.end = 0, self.ndata-1
      self.times = numpy.arange(self.ndata)
      self.calculate()
    def calculate(self):
      self.av = self.values[self.start:self.end+1].mean()
      self.var = self.values[self.start:self.end+1].var()
      self.sigma = math.sqrt(self.var)
      #Compute autocorrelation
      aclen = (self.end-self.start+1)/4
      if aclen>1:
        a = self.values[self.start:self.end+1] - self.av
        a = numpy.correlate(a,a[aclen:-aclen],'valid')
        self.autocorr = (a[aclen:]+a[aclen::-1])/(2.*a[aclen])
        cut = numpy.where(self.autocorr < 0.)[0]
        if len(cut)>0:
          cut = cut[0]
          if 5*cut<len(self.autocorr): self.autocorr = self.autocorr[:5*cut]
        self.actime = 2*sum(self.autocorr)-1
        if self.actime < 1.: self.actime=1.
      else:
        self.autocorr = numpy.array([1.,0.])
        self.actime = 1.
      self.err = math.sqrt(self.var*self.actime/(self.end-self.start))
      #Do blocking analysis.
      nblock = int(math.log(self.end-self.start+1,2)-3)
      self.blocking = numpy.zeros(nblock+1)
      data = self.values[self.start:self.end+1]
      for i in range(0,nblock+1):
        self.blocking[i] = math.sqrt(data.var()/(len(data)-1.))
        data = 0.5*(data[:-1:2]+data[1::2])

  class TracePlotWidget(MyMplCanvas):
    def __init__(self, model):
      self.model = model
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.20,0.15,0.67,0.82])
      self.histAxes = self.figure.add_axes([0.87,0.15,0.12,0.82])
      self.plotFigure()
    def plotFigure(self):
      ax = self.axes; hax = self.histAxes
      ax.cla(); hax.cla()
      plt.rc('font', family='serif', size=9)
      # Plot trace.
      times,values = self.model.times, self.model.values
      start,end = self.model.start, self.model.end
      ax.plot(times[start:end+1],values[start:end+1],color='b')
      av,err = self.model.av, self.model.err
      ax.axhspan(av-err, av+err, color='y', alpha=0.4)
      ax.axhline(y=av, color='k', lw=0.2)
      ax.axis(xmin=start, xmax=end)
      ax.set_xlabel(r"Monte Carlo timestep")
      ax.set_ylabel(self.model.scalar.name)
      # Plot histogram.
      hax.hist(self.model.values[start:end+1],
               bins=32, orientation='horizontal',
               range=ax.get_ybound())
      hax.set_ybound(ax.get_ybound())
      hax.set_xticks([]); hax.set_yticks([])
       
  class AutocorrPlotWidget(MyMplCanvas):
    def __init__(self, model):
      self.model = model
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      ax = self.axes = self.figure.add_axes([0.15,0.15,0.82,0.82])
      self.plotFigure()
    def plotFigure(self):
      plt.rc('font', family='serif', size=9)
      ax = self.axes
      ax.cla()
      autocorr = self.model.autocorr
      ax.axhline(y=0, color='k', lw=0.3)
      ax.plot(numpy.arange(len(autocorr)),autocorr)
      ax.axis(xmax=len(autocorr)-1)
      ax.set_xlabel(r"simulation time")
      ax.set_ylabel(r"autocorrelation")

  class BlockingPlotWidget(MyMplCanvas):
    def __init__(self, model):
      self.model = model
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.15,0.15,0.82,0.82])
      self.plotFigure()
    def plotFigure(self):
      plt.rc('font', family='serif', size=9)
      ax = self.axes
      ax.cla()
      blocking = self.model.blocking
      ax.plot(numpy.arange(len(blocking)),blocking)
      ax.axhline(y=self.model.err, color='k', lw=1.0)
      ax.axis(ymin=0)
      ax.set_xlabel(r"blocking step")
      ax.set_ylabel(r"error estimate")
