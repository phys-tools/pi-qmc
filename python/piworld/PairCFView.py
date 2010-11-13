#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *


class PairCFView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    pylab.rc('font', family='serif', size=9)

    self.pcf = data.file.getPairCF(estimatorNode.name)#,Unit.nm)

    self.plot =  self.Plot1DWidget(self)
    vbox = QtGui.QVBoxLayout()
    vbox.setSpacing(0)
    vbox.addWidget(self.plot,10)

    self.setLayout(vbox)
    self.setEnabled(True)

  class Plot1DWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      self.pcf = data.pcf.data
      self.error = data.pcf.error
      self.r, self.dr = data.pcf.getCoordinateGrid(4.7595*52.9177)
      self.vol = 4./3.*math.pi*(
                  (self.r+0.5*self.dr)**3-(self.r-0.5*self.dr)**3)
      self.rmax = self.r[-1]+0.5*self.dr
      #self.origin = self.data.density.origin
      #self.scale = self.data.density.scale
      #self.shape = numpy.shape(self.rho)
      #self.limit = numpy.array([0,self.shape[0],0,self.shape[1]]) 
      #self.extent =  numpy.array(
      #    [self.origin[0]+0.5*self.scale[0],
      #     self.origin[0]+(self.shape[0]-0.5)*self.scale[0],
      #     self.origin[1]+0.5*self.scale[1],
      #     self.origin[1]+(self.shape[1]-0.5)*self.scale[1]])
      #self.currentExtent = self.extent.copy()
      #aspect = (self.extent[3]-self.extent[2])/(self.extent[1]-self.extent[0]) 
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.01,0.15,0.98,0.83])
      self.plotFigure()
    def plotFigure(self):
      self.axes.cla()
      #self.axes.imshow(self.rho[self.limit[0]:self.limit[1],
      #                          self.limit[2]:self.limit[3]].transpose(),
      #    origin='lower', interpolation='bicubic', extent=self.currentExtent)
      
      self.axes.plot(self.r,self.pcf/self.vol)
      if self.error != None:
        self.axes.errorbar(self.r, self.pcf/self.vol,
                                   self.error/self.vol)
      self.axes.axis(xmax=self.rmax)
      self.axes.set_yticks([])
      self.axes.set_xlabel(r"$x$ (pm)")

