"""Create object tree from pi-qmc files


Classes:

    File

Functions:

    openFile([name])

"""

import tables
from units import Unit
from scalar import Scalar
from response import ResponseFunction

def openFile(name="pimc.h5"):
  return File(name)

class File(object):
  """
  In-memory representation of a pimc.h5 file.
  """

  def __init__(self,name):
    self.file=tables.openFile(name)

  def close(self):
    self.file.close()

  def __str__(self):
    return self.file.__str__()


  def getTemperature(self,unit=Unit.Ha):
    return unit.out(self.file.getNode("/simInfo","temperature").read()[0])

  def getNSlice(self):
    return int(self.file.getNode("/simInfo","nslice").read()[0])

  def getNPart(self):
    return int(self.file.getNode("/simInfo","npart").read()[0])

  def getSuperCell(self,unit=Unit.a0):
    return unit.out(self.file.getNode("/simInfo","superCell").read())

  def getScalar(self,name,unit=None):
    node = self.file.getNode("/estimators",name)
    try:
      unitName = node.getAttr("unit")
      dataUnit = Unit.getByName(unitName[0])
    except:
      if unit == None:
        dataUnit = None
      else:
        dataUnit = unit.native
    nstep = self.file.getNodeAttr("/estimators","nstep")[0]
    if unit == None:
      data = node.read()[:nstep]
    else:
      data = unit.convert(node.read()[:nstep],dataUnit)
    scalar = Scalar(name,data,unit)
    return scalar

  def getResponseFunction(self,name):
    data = self.file.getNode("/estimators",name).read()
    error = self.file.getNode("/estimators",name+"_err").read()
    t = self.file.getNode("/simInfo","temperature").read()[0]
    return ResponseFunction(name,data,error,t)

