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
from density import Density
from paircf import PairCF

class estimatorNode(object):
  """
  Representation of a node in the estimator group of a pimc.h5 file.
  """
  def __init__(self,node,typeString=""):
    self.name = node.name
    self.typeString = typeString
    self.node = node
    # Hack to add type strings to old files.
    if typeString=="":
      if self.name[:4]=="perm":
        self.typeString="histogram/permutation"
      if self.name[:3]=="rho":
        if self.name[-3:]=="(r)":
          self.typeString="array/pair-correlation"
        else:
          self.typeString="array/density"
      if self.name[:12]=="conductivity":
        self.typeString="dynamic-array/conductivity"
      if self.name[:11]=="conductance":
        self.typeString="dynamic-array/conductance"
      if self.name[-6:]=="energy":
        self.typeString="scalar-energy"
      if self.name=="coulomb_energy":
        self.typeString="scalar-energy/coulomb-energy"
      if self.name=="thermo_energy":
        self.typeString="scalar-energy/thermo-energy"
      if self.name=="virial_energy":
        self.typeString="scalar-energy/virial-energy"
      if self.name=="free_energy":
        self.typeString="histogram/free-energy"
      if self.name[:11]=="bond_length":
        self.typeString="scalar-length/bond-length"
      if self.typeString=="" and (self.name[:1]=="g"or 
                                  self.name=="bondlength"):
        self.typeString="array/pair-correlation"
      if self.typeString=="" and self.name[:12]=="bond_length_":
        self.typeString="scalar-length/bond-length"
      if self.name=="winding":
        self.typeString="histogram/winding"
      if self.name=="chiM":
        self.typeString="scalar/diamagnetism"
      if self.name[:13]=="dipole_moment":
        self.typeString="scalar"
    if self.name[-4:]=="_err":
      self.typeString += "-err" 
  

def openFile(name="pimc.h5"):
  return File(name)

class File(object):
  """
  Representation of a pimc.h5 file.
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
    try:
      error = self.file.getNode("/estimators",name+"_err").read()
    except:
      error = None
    t = self.file.getNode("/simInfo","temperature").read()[0]
    return ResponseFunction(name,data,error,t)

  def getDensity(self,name,unit=Unit.a0):
    node = self.file.getNode("/estimators",name)
    data = node.read()
    try:
      error = self.file.getNode("/estimators",name+"_err").read()
    except:
      error = None
    origin = unit.convert(node.getAttr("origin"),unit.a0)
    scale = unit.convert(node.getAttr("scale"),unit.a0)
    return Density(name,data,error,origin,scale)

  def getPairCF(self,name):
    node = self.file.getNode("/estimators",name)
    data = node.read()
    try:
      error = self.file.getNode("/estimators",name+"_err").read()
    except:
      error = None
    return PairCF(name,data,error)

  def getAllSpecies(self):
    slist = []
    return slist

  def getEstimatorList(self):
    elist = []
    try:
      for node in self.file.iterNodes("/estimators"):
        elist.append(estimatorNode(node))
    except tables.NoSuchNodeError:
      pass
    return elist
