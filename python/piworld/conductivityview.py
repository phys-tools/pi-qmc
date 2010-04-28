#-*- coding: utf-8 -*-
import math,scipy,numpy,tables,pylab
from scipy.fftpack import *
from scipy.optimize import leastsq
import tables, pitools
from pitools import Unit
from estimatorview import *


class ConductivityView(EstimatorView):

  def __init__(self, estimatorNode, data, parent=None):
    EstimatorView.__init__(self,parent)

    h5file = data.file.file
    jj = h5file.getNode("/estimators","conductivity").read()
    jj_err = h5file.getNode("/estimators","conductivity_err").read()
    #ej = h5file.getNode("/estimators","eind").read()
    #ej_err = h5file.getNode("/estimators","eind_err").read()
    #l = h5file.getNode("/simInfo","superCell").read()[0]*0.05291772108
    l = data.superCell[0]
    npart = data.file.getNPart()
    nslice = data.nslice
    temperature = data.file.getTemperature(Unit.Ha)

    self.ndx,self.nx,self.nfreq = jj.shape
    deltax=l/0.05291771208/self.nx
    density=npart/(l/0.05291771208)

    #print "  Correlation function array has dimensions " + \
    #    ("nx=%i, ndx=%i, nfreq=%i." % (self.nx,self.ndx,self.nfreq))

    #Calculate the conductance.
    jjw=numpy.mean(jj,1)
    jjw_err=numpy.mean(jj_err,1)
    self.omega1=2.*numpy.pi*temperature
    self.omegan=self.omega1*numpy.array(range(1,self.nfreq))
    self.sigma0iw = jjw[:,1:]*2.*numpy.pi/self.omegan
    self.sigma0iw = scipy.fftpack.fftshift(self.sigma0iw,[0])
    self.sigma0iw_err = jjw_err[:,1:]*2.*numpy.pi/self.omegan
    self.sigma0iw_err = scipy.fftpack.fftshift(self.sigma0iw_err,[0])

    #Fit sigma
    p0=[0.86,4.*self.omega1,1.0,20.*self.omega1]
    plsq = leastsq(self.residuals,p0,
                   args=(self.sigma0iw[self.ndx/2,1:],
                         self.sigma0iw_err[self.ndx/2,1:],
                         self.omegan[1:]),maxfev=2000)
    self.pfit_a=plsq[0]
    #print "\n Fitting parameters: σ_0=%f, ℏω_0=%f eV, ℏω_0=%f eV"\
    #      %(self.peval(0,self.pfit_a),
    #        self.pfit_a[1]*27.211,self.pfit_a[3]*27.211)


    pylab.rc('font', family='serif', size=9)

    self.plot =  self.PlotWidget(self)
    vbox = QtGui.QVBoxLayout()
    vbox.setSpacing(0)
    vbox.addWidget(self.plot,10)

    self.info = self.InfoWidget(self)
    vbox.addWidget(self.info,1)

    self.setLayout(vbox)
    


#    pylab.clf()
#    pylab.axes([0.15,0.15,0.83,0.83])
#    #pylab.plot(omegan/omega1,sigmaiw[ndx/2,:nfreq])
#    pylab.errorbar(omegan[:40]/omega1,sigma0iw[ndx/2,:40],
#                   sigma0iw_err[ndx/2,:40], fmt=".")
#    pylab.plot(omegas/omega1,peval(omegas,pfit_a))
#    pylab.axis([0,20.5,0,2.5])
#    pylab.xlabel(r"$\omega_n/\omega_1$")
#    pylab.ylabel(r"$\sigma^0(x,x;i\omega_n)\ (e^2/h)$")
#    pylab.savefig("sigma0iwn.pdf")
#    pylab.savefig("sigma0iwn.png",dpi=300)
    
#pylab.clf()
##pylab.plot(omegan/omega1,sigmaiw[ndx/2,:nfreq]/viw)
#pylab.errorbar(omegan/omega1,sigmaiw[ndx/2,:],sigmaiw_err[ndx/2,:], fmt=".")
#pylab.axis([0,nfreq,0,2.5])
#pylab.xlabel(r"$\omega_n/\omega_1$")
#pylab.ylabel(r"$\sigma(x,x;i\omega_n)\ (G_0)$")
#pylab.savefig("sigmaiwn.pdf")
##pylab.savefig("sigmaiwn.png")

#pylab.clf()
#pylab.hold(True)
#pylab.errorbar(xarray,vindiw[:,0],vindiw_err[:,0])
#pylab.plot(xarray0,vappiw[:,0])
#pylab.errorbar(xarray0,vtotiw[:,0],vtotiw_err[:,0])
#pylab.xlabel(r"$x\ (\rm{nm})$")
#pylab.savefig("vindiwn.pdf")
#pylab.savefig("vindiwn.png")
#
#    pylab.clf()
#    pylab.axes([0.18,0.15,0.75,0.72])
#    pylab.imshow(numpy.log(numpy.abs(sigma0iw[ndx/4:-ndx/4,:21])),origin="lower",
#                 extent=[1,21,-l/4,l/4], aspect='auto',
#                 interpolation="spline16",vmin=-8,vmax=math.log(sigma0iw.max()))
#    pylab.xlabel(r"$\omega_n/\omega_1$")
#    pylab.ylabel(r"$\Delta x\ \rm{(nm)}$")
#    pylab.title(r"$-\chi^0_{j_x j_x}(x+\Delta x,x;i\omega_n)/\omega_n$")
#    pylab.jet()
#    pylab.colorbar()
#    pylab.savefig("jj0dxiwn.pdf")
#    pylab.savefig("jj0dxiwn.png",dpi=300)
    
#pylab.clf()
#pylab.hold(True)
#pylab.imshow(sigmaiw/viw,origin="lower",extent=[1,nfreq,-l/2,l/2],
#             aspect=nfreq/l, interpolation="spline16",vmin=0,vmax=1.5)
#pylab.xlabel(r"$\omega_n/\omega_1$")
#pylab.ylabel(r"$\Delta x\ \rm{(nm)}$")
#pylab.title(r"$-\chi_{j_x j_x}(x+\Delta x,x;i\omega_n)/\omega_n$",font)
#pylab.jet()
#pylab.colorbar()
#pylab.savefig("jjdxiwn.pdf")
#pylab.savefig("jjdxiwn.png")
#pylab.hold(False)
#
#file = open("giwn.dat","w")
#for i in range(nfreq-1):
#  file.write("%10.6f %10.6f\n"%(omegan[i],sigmaiw[ndx/2,i]))

  #Fit to free-particle response function.
  def peval(self,omegan,p):
    #return p[0]/(omegan**2+p[1]**2)+p[2]/(omegan**2+p[3]**2)
    return p[0]/numpy.sqrt(0.5*(1+numpy.sqrt(1+omegan**2/p[1]**2)))\
          +p[2]/numpy.sqrt(0.5*(1+numpy.sqrt(1+omegan**2/p[3]**2)))
    #      +p[2]/(1+omegan**2/p[3]**2)
  def residuals(self,p,sigmaiwn,sigmaiwn_err,omegan):
    return (sigmaiwn-self.peval(omegan,p))/sigmaiwn_err

  class PlotWidget(MyMplCanvas):
    def __init__(self, data):
      self.data = data
      MyMplCanvas.__init__(self)
    def computeInitialFigure(self):
      maxfreq = math.ceil((self.data.omegan/self.data.omega1).max())
      maxg = math.ceil((self.data.sigma0iw).max())
      omegas = numpy.arange(0,maxfreq+0.4,0.1)*self.data.omega1
      self.axes = self.figure.add_axes([0.15,0.15,0.83,0.83])
      #pylab.plot(omegan/omega1,sigmaiw[ndx/2,:nfreq])
      self.axes.errorbar(self.data.omegan[:]/self.data.omega1, 
                         self.data.sigma0iw[self.data.ndx/2,:],
                         self.data.sigma0iw_err[self.data.ndx/2,:], fmt=".")
      self.axes.plot(omegas/self.data.omega1,
                     self.data.peval(omegas,self.data.pfit_a))
      self.axes.axis([0,maxfreq+0.4,0,maxg+0.4])
      self.axes.set_xlabel(r"$\omega_n/\omega_1$")
      self.axes.set_ylabel(r"$G$ ($e^2/h$)")

  class InfoWidget(QtGui.QWidget): 
    def __init__(self, data, parent=None):
      self.data = data
      QtGui.QWidget.__init__(self,parent)
 
      gfit = self.data.peval(0.,self.data.pfit_a)

      hbox = QtGui.QHBoxLayout()
      hbox.setSpacing(3)
      self.labelGFit = QtGui.QLabel(u"G = %g ± ? (e<sup>2</sup>/h)" % gfit)
      hbox.addWidget(self.labelGFit)
      hbox.addStretch(1)
      self.setLayout(hbox)
