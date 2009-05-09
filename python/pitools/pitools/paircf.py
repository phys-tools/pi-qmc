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


