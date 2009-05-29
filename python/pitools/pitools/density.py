"""Density storage and analysis


Classes:

    Density


"""
import numpy,math

class Density(object):

  def __init__(self,name,data,error,origin,scale):
    self.name = name
    self.data = data
    self.error = error
    self.origin = origin
    self.scale = scale
    self.extent = scale*data.shape

  def __str__(self):
    return "Density '%s'" % (self.name)


