#!/usr/bin/env python

from subprocess import Popen, PIPE, STDOUT, CalledProcessError, call
import xml.etree.ElementTree
import os, sys
import tables

def runHarness():
    greeting()
    check("ThermalEnergyEstimator", "thermo_energy")
    check("SpinEnergyEstimator", "")
    check("VirialEnergyEstimator", "virial_energy")
    check("CoulombEnergyEstimator", "coulomb_energy")
    check("AngularMomentumEstimator", "angular_momentum")
    check("BondLengthEstimator", "bond_length_e-e")
    check("FrequencyEstimator", "frequency")
    check("BoxEstimator", "box_z_e")
    check("PositionEstimator", "position_z_e")
    check("ConductivityEstimator", "conductivity")
    check("ConductivityEstimator2D", "conductivity2D")
    check("ConductanceEstimator", "conductance")
    check("DensityEstimator", "rhoe")
    check("DensDensEstimator", "nne")
    check("DensCountEstimator", "counte")
    check("CountCountEstimator", "counte")
    check("SpinChoiceDensityEstimator", "rhoeup")
    check("PairCFEstimator", "geh")
    check("ZeroVarDensityEstimator", "density")
    check("SpinChoicePCFEstimator", "geh")
    check("PermutationEstimator", "perm_e")
    check("JEstimator", "exchange")
    check("DiamagneticEstimator", "chiM")
    check("WindingEstimator", "winding")
    check("SKOmegaEstimator", "skomega")
    check("EMARateEstimator", "ema_rate")

def greeting():
    print "\nPython test harness for pi\n"
    print "  Note: This script runs a bunch of input files for pi and verifies"
    print "        that the EstimatorParser continues to work " + \
        "during refactoring."
    print

def check(tagName, outputName):
    print "Testing <%s> estimator tag" % tagName
    createXMLInputFile(tagName);
    piexe = sys.argv[1]
    try:
        command  = [piexe, "%s.xml" % tagName]
        process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if exitsWithError(process):
            stderr = process.stderr.read()
            print stderr,
            reportFailed("process failed <%s>" % tagName)
            sys.exit()
        else:
            checkPimcH5Output(outputName) 
    except tables.exceptions.NoSuchNodeError as details:
        print details
        os.system("h5ls pimc.h5/estimators")
        reportFailed("<%s>" % tagName)
        sys.exit()
    else:
        reportPassed("<%s>" % tagName)

def checkPimcH5Output(outputName):
    if (outputName != ""):
        h5file = tables.openFile("pimc.h5")
        node = h5file.getNode("/estimators", outputName)
        node.read()
        h5file.close()

def exitsWithError(process):
  return process.wait() != 0

RED = "\033[91m"
GREEN = "\033[92m"
ENDCOLOR = "\033[0m"

def redText(string):
    return RED + string + ENDCOLOR

def greenText(string):
    return GREEN + string + ENDCOLOR

def reportPassed(string=""):
    print greenText("...PASSED  " + string)

def reportFailed(string=""):
    print redText("...FAILED  " + string)

def createXMLInputFile(tagName):
    xmlString = \
        '<?xml version="1.0"?>' + \
        '<Simulation>' + \
        '  <SuperCell a="25 A" x="1" y="1" z="1"/>' + \
        '  <Species name="e" count="1" mass="1 m_e" charge="-1"/>' + \
        '  <Species name="h" count="1" mass="1 m_e" charge="-1"/>' + \
        '  <Temperature value="1 Ha" tau="0.1 Ha-1"/>' + \
        '  <Action>' + \
        '    <SpinChoiceFixedNodeAction/>' + \
        '  </Action>' + \
        '  <Estimators/>' + \
        '  <PIMC>' + \
        '    <RandomGenerator/>' + \
        '    <Measure estimator="all"/>' + \
        '    <Collect estimator="all"/>' + \
        '  </PIMC>' + \
        '</Simulation>'
    xmlTree = xml.etree.ElementTree.ElementTree(
                  xml.etree.ElementTree.fromstring(xmlString));
    estimatorsNode = xmlTree.find("Estimators")
    estimatorNode = xml.etree.ElementTree.Element(tagName)
    addSpecialAttributes(estimatorNode)
    estimatorsNode.append(estimatorNode)
    fileName = tagName + ".xml"
    xmlTree.write(fileName, "utf-8", True)
    os.system("xmllint --format %s > pimc.xml" % fileName)
    os.system("cp pimc.xml %s" % fileName)

def addSpecialAttributes(node):
    if (node.tag == "DensityEstimator" or \
            node.tag == "DensDensEstimator" or \
            node.tag == "DensCountEstimator" or \
            node.tag == "CountCountEstimator" or \
            node.tag == "SpinChoiceDensityEstimator"):
        node.set("nx", "4")
        node.set("ny", "4")
        node.set("nz", "4")
        node.set("a", "1.0")
        node.set("species", "e")
    elif (node.tag == "JEstimator"):
        node.set("nBField", "1")
    elif (node.tag == "SKOmegaEstimator"):
        node.set("nx", "4")
        node.set("ny", "4")
        node.set("nz", "4")
        node.set("nfreq", "1")
    elif (node.tag == "PairCFEstimator" or
          node.tag == "SpinChoicePCFEstimator"):
        node.set("species1", "e")
        node.set("species2", "h")
        node.set("name", "geh")
        radialNode = xml.etree.ElementTree.Element("Radial")
        radialNode.set("nbin", "10")
        radialNode.set("max", "12.5 A")
        node.append(radialNode)
    elif (node.tag == "ZeroVarDensityEstimator"):
        node.set("nspecies", "1")
        node.set("name", "density")
       



runHarness()
