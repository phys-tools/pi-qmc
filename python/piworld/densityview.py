#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from estimatorview import *


class DensityView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    self.density = data.file.getDensity(estimatorNode.name,Unit.nm)

    if self.density.rank == 2:
      self.plot =  self.Plot2DWidget(self)
    vbox = QtGui.QVBoxLayout()
    vbox.setSpacing(0)
    vbox.addWidget(self.plot,10)

    self.setLayout(vbox)
    
  class Plot2DWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      self.rho = data.density.data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.18,0.15,0.80,0.83])
      origin = self.data.density.origin
      extent = self.data.density.extent
      scale = self.data.density.scale
      shape = numpy.shape(self.rho)
      self.axes.imshow(self.rho.transpose(), origin='lower',
          interpolation='bicubic',
          extent=(origin[0]+0.5*scale[0], origin[0]+(shape[0]-0.5)*scale[0],
                  origin[1]+0.5*scale[1], origin[1]+(shape[1]-0.5)*scale[1]) )
      self.axes.set_xlabel(r"$x$ (nm)")
      self.axes.set_ylabel(r"$y$ (nm)")

  class Plot3DWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.18,0.15,0.80,0.83])
      origin = self.data.density.origin
      extent = self.data.density.extent
      scale = self.data.density.scale
      shape = numpy.shape(rho)
      self.axes.imshow(self.rho.sum(3).transpose(), origin='lower',
          interpolation='bicubic',
          extent=(origin[0]+0.5*extent[0], origin[1]+0.5*extent[1],
                  origin[0]+(shape[0]-0.5)*extent[0],
                  origin[1]+(shape[1]-0.5)*extent[1]) )
