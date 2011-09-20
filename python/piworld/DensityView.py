#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from EstimatorView import *


class DensityView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    pylab.rc('font', family='serif', size=9)

    self.grabGesture(QtCore.Qt.PanGesture)
    self.grabGesture(QtCore.Qt.PinchGesture)
    self.grabGesture(QtCore.Qt.SwipeGesture)

    self.density = data.file.getDensity(estimatorNode.name,Unit.nm)

    if self.density.rank == 2:
      self.plot =  self.Plot2DWidget(self)
    elif self.density.rank == 3:
      self.plot =  self.Plot3DWidget(self)
    vbox = QtGui.QVBoxLayout()
    vbox.setSpacing(0)
    vbox.addWidget(self.plot,10)

    self.setLayout(vbox)
    self.setEnabled(True)

  def mousePressEvent(self, event):
    print "Mouse press event"

  def wheelEvent(self, event):
    print "Wheel event"

  def event(self, event):
    if event.type() == QtCore.QEvent.Gesture:
      pan = event.gesture(QtCore.Qt.PanGesture)
      pinch = event.gesture(QtCore.Qt.PinchGesture)
      swipe = event.gesture(QtCore.Qt.SwipeGesture)
      if swipe: self.plot.swipeTriggered(swipe)
      return True
    return False
    
  class Plot2DWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      self.rho = data.density.data
      self.origin = self.data.density.origin
      self.scale = self.data.density.scale
      self.shape = numpy.shape(self.rho)
      self.limit = numpy.array([0,self.shape[0],0,self.shape[1]]) 
      self.extent =  numpy.array(
          [self.origin[0]+0.5*self.scale[0],
           self.origin[0]+(self.shape[0]-0.5)*self.scale[0],
           self.origin[1]+0.5*self.scale[1],
           self.origin[1]+(self.shape[1]-0.5)*self.scale[1]])
      self.currentExtent = self.extent.copy()
      aspect = (self.extent[3]-self.extent[2])/(self.extent[1]-self.extent[0]) 
      if aspect > 3: #zoom in y direction
        jmid = self.shape[1]/2
        self.limit[3] = jmid + self.shape[0]
        self.limit[2] = jmid - self.shape[0]
        self.currentExtent[3] = self.origin[1] \
                              +(self.limit[3]-0.5)*self.scale[1]
        self.currentExtent[2] = self.origin[1] \
                              +(self.limit[2]-0.5)*self.scale[1]
      if aspect < 0.3: #zoom in x direction
        imid = self.shape[0]/2
        self.limit[1] = imid + self.shape[1]
        self.limit[0] = imid - self.shape[1]
        self.currentExtent[1] = self.origin[0] \
                              +(self.limit[1]-0.5)*self.scale[0]
        self.currentExtent[0] = self.origin[0] \
                              +(self.limit[0]-0.5)*self.scale[0]
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.18,0.15,0.80,0.83])
      self.plotFigure()
    def plotFigure(self):
      self.axes.cla()
      self.axes.imshow(self.rho[self.limit[0]:self.limit[1],
                                self.limit[2]:self.limit[3]].transpose(),
          vmin = 0.,
          origin='lower', interpolation='bicubic', extent=self.currentExtent)
      # Show lighter bars at edges when zoomed.
      if self.limit[0]>0:
        self.axes.axvspan(
             self.currentExtent[0],
             0.980*self.currentExtent[0]+0.020*self.currentExtent[1],
             color='w',alpha=0.2)
      if self.limit[1]<self.shape[0]-1:
        self.axes.axvspan(
             0.020*self.currentExtent[0]+0.980*self.currentExtent[1],
             self.currentExtent[1], color='w',alpha=0.2)
      if self.limit[2]>0:
        self.axes.axhspan(self.currentExtent[2],
             0.980*self.currentExtent[2]+0.020*self.currentExtent[3],
             color='w',alpha=0.2)
      if self.limit[3]<self.shape[1]-1:
        self.axes.axhspan(
             0.020*self.currentExtent[2]+0.980*self.currentExtent[3],
             self.currentExtent[3], color='w',alpha=0.2)
      self.axes.axis(self.currentExtent)
      self.axes.set_xlabel(r"$x$ (nm)")
      self.axes.set_ylabel(r"$y$ (nm)")

    def swipeTriggered(self,swipe):
      hdir = swipe.horizontalDirection()
      vdir = swipe.verticalDirection()
      if hdir == swipe.Right and self.limit[0]>0:
        imin = (12*self.limit[0]-2*self.limit[1])/10
        if imin < 0: imin = 0
        self.limit[1] += imin-self.limit[0]
        self.limit[0] = imin
        self.currentExtent[1] = self.origin[0] \
                              +(self.limit[1]-0.5)*self.scale[0]
        self.currentExtent[0] = self.origin[0] \
                              +(self.limit[0]-0.5)*self.scale[0]
        self.plotFigure()
        self.draw()
      if hdir == swipe.Left and self.limit[1]<self.shape[0]-1:
        imax = (12*self.limit[1]-2*self.limit[0])/10
        if imax > self.shape[0]-1: imax = self.shape[0]-1
        self.limit[0] += imax-self.limit[1]
        self.limit[1] = imax
        self.currentExtent[1] = self.origin[0] \
                              +(self.limit[1]-0.5)*self.scale[0]
        self.currentExtent[0] = self.origin[0] \
                              +(self.limit[0]-0.5)*self.scale[0]
        self.plotFigure()
        self.draw()

  class Plot3DWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      self.rho = data.density.data
      self.origin = self.data.density.origin
      self.scale = self.data.density.scale
      self.shape = numpy.shape(self.rho)
      self.limit = numpy.array([0,self.shape[0],0,self.shape[1]]) 
      self.extent =  numpy.array(
          [self.origin[0]+0.5*self.scale[0],
           self.origin[0]+(self.shape[0]-0.5)*self.scale[0],
           self.origin[1]+0.5*self.scale[1],
           self.origin[1]+(self.shape[1]-0.5)*self.scale[1]])
      self.currentExtent = self.extent.copy()
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      self.axes = self.figure.add_axes([0.18,0.15,0.80,0.83])
      self.plotFigure()
    def plotFigure(self):
      self.axes.imshow(self.rho.sum(2).transpose(), origin='lower',
          interpolation='bicubic', extent=self.currentExtent)
      self.axes.axis(self.currentExtent)
      self.axes.set_xlabel(r"$x$ (nm)")
      self.axes.set_ylabel(r"$y$ (nm)")

    def swipeTriggered(self,swipe):
      pass
