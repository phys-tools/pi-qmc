"""Response function storage and analysis


Classes:

    ResponseFunction


"""
import numpy,math

class ResponseFunction(object):

  def __init__(self,name,data,error,temperature):
    self.name = name
    self.data = data
    self.error = error
    self.nfreq = data.shape[-1]
    self.omega1 = 2*math.pi*temperature
    self.omega = numpy.arange(0,self.nfreq,1)*self.omega1

  def __str__(self):
    return "Response function '%s'" % (self.name)


