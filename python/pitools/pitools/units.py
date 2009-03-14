"""Unit conversion utilities.


Classes:

    Unit


"""

class Unit(object):
  allUnits = {}
  def __init__(self,name=""):
    self.name = name
    self.conversion = 1
    self.native = self
  def __str__(self):
    return self.name
  def out(self,x):
    return self.conversion*x
  def convert(self,x,fromUnit):
    return self.conversion/fromUnit.conversion*x
  def getByName(name):
    return Unit.allUnits[name]
  getByName = staticmethod(getByName)
  def addToDictionary(unit):
    Unit.allUnits[unit.name] = unit
  addToDictionary = staticmethod(addToDictionary)

class EnergyUnit(Unit):
  Ha = Unit("Ha")
  def __init__(self,name,conversion):
    self.name = name
    self.conversion = conversion
    self.native = EnergyUnit.Ha
    Unit.addToDictionary(self)

class LengthUnit(Unit):
  a0 = Unit("a0")
  def __init__(self,name,conversion):
    self.name = name
    self.conversion = conversion
    self.native = LengthUnit.a0
    Unit.addToDictionary(self)

Unit.Ha = EnergyUnit.Ha
Unit.K = EnergyUnit("K",3.1577465e5)     
Unit.mK = EnergyUnit("mK",3.1577465e8)
Unit.uK = EnergyUnit("uK",3.1577465e11)
Unit.nK = EnergyUnit("nK",3.1577465e14)
Unit.eV = EnergyUnit("eV",27.211396)
Unit.meV = EnergyUnit("meV",27211.396)

Unit.a0 = LengthUnit.a0
Unit.pm = LengthUnit("pm",52.9177)
Unit.A = LengthUnit("A",0.529177)
Unit.nm = LengthUnit("nm",0.0529177)
Unit.um = LengthUnit("um",5.29177e-5)
