"""Scalar dataset storage and analysis


Classes:

    Scalar


"""
import numpy,math

class Scalar(object):

  def __init__(self,name,data,unit):
    self.name = name
    self.data = data
    self.unit = unit

  def __str__(self):
    return "Scalar '%s' (%s)" % (self.name,self.unit.name)

  def getAverage(self):
    av = self.data.mean()
    npts = self.data.size
    var = self.data.var()
    err = math.sqrt(var/(npts-1))
    return av,err

