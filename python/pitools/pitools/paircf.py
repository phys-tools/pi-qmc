"""Pair correlation function storage and analysis


Classes:

    PairCF


"""
import numpy,math

class PairCF(object):

  def __init__(self,name,data,error):
    self.name = name
    self.data = data
    self.error = error

  def __str__(self):
    return "Pair correlation function '%s'" % (self.name)

  def getCoordinateGrid(self, rmax, idim=0, rmin=0.):
    ngrid = self.data.shape[idim]
    dr = (rmax-rmin)/ngrid
    r = (numpy.arange(ngrid)+0.5)*dr+rmin
    return r,dr;
