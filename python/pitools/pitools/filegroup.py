"""Create parallel object trees from pi-qmc files


Classes:

    FileGroup

Functions:

    openFile([name])

"""

import tables, numpy, glob
from file import File
from units import Unit
from scalar import Scalar
from response import ResponseFunction
from density import Density
from paircf import PairCF

def openFileGroup(nameList=None,pattern=None):
  if pattern!=None:
    if nameList:
      nameList.append(glob.glob(pattern))
    nameList = glob.glob(pattern)
  fileList = []
  for name in nameList:
    fileList.append(File(name))
  return FileGroup(fileList)

class FileGroup(object):
  """
  Representation of a pimc.h5 file.
  """

  def __init__(self,fileList):
    self.fileList = fileList
    self.nfile = len(fileList)

  def close(self):
    for file in self.fileList:
      file.close()

  def __str__(self):
    return "Filegroup"

  def getNSlice(self):
    nslice = numpy.zeros(self.nfile,numpy.int32)
    for i,file in enumerate(self.fileList):
      nslice[i] = file.getNSlice()
    return nslice

  def getScalar(self,name,unit=None):
    scalar = numpy.zeros(self.nfile)
    error = numpy.zeros(self.nfile)
    for i,file in enumerate(self.fileList):
      scalar[i],error[i] = file.getScalar(name,unit).getAverage()
    return scalar,error

