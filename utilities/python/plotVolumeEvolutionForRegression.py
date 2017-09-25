import random as rd
import vtk
from vtk import vtkPolyDataReader, vtkBooleanOperationPolyDataFilter, vtkMassProperties, vtkPolyDataWriter, vtkIntersectionPolyDataFilter
import numpy as np
import sys
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import os


## Warning : the computed volumes are somewhat unstable when meshes are highly non-convex.

def extractIndex(filename):
    for k2 in range(1, len(filename) + 1):
        char = filename[- k2]
        if char == '.':
            k2 = k2 - 1
            break
    for k1 in range(k2 + 1, len(filename) - k2):
        char = filename[- k1]
        if char == '_':
            k1 = k1 - 1
            break
    index = int(float(filename[-k1:-k2]))
    return index

def extractInfoFromModel(pathToRegressionFolder):
    structureNames = []
    modelTree = ET.parse(os.path.join(pathToRegressionFolder, "model.xml")).getroot()
    deformationParameters = modelTree.find('deformation-parameters')
    t0 = float(deformationParameters.find('t0').text)
    tn = float(deformationParameters.find('tn').text)
    nbTimePoints = int(deformationParameters.find('number-of-timepoints').text)
    template = modelTree.find('template')
    objects = template.findall('object')
    for obj in objects:
        structureNames.append(obj.attrib['id'])
    structureNames.sort()
    return structureNames, t0, tn, nbTimePoints

def extractInfoFromDataSet(pathToRegressionFolder, structureNames):
    dataMeshPaths = []
    dataAges = []
    dataSetTree = ET.parse(os.path.join(pathToRegressionFolder, "data_set.xml")).getroot()
    subjects = dataSetTree.findall('subject')
    if len(subjects)>1: # OLD XML CONVENTION : TO BE REMOVED SOMEDAY SOON
        for subject in subjects:
            visit = subject.find('visit')
            dataAges.append(float(visit.find('age').text))
        for struct in structureNames:
            structDataMeshPaths = []
            for subject in subjects:
                visit = subject.find('visit')
                filenames = visit.findall('filename')
                for filename in filenames:
                    objectId = filename.attrib['object_id']
                    if objectId == struct:
                        structDataMeshPaths.append(os.path.normpath(os.path.join(pathToRegressionFolder, "output", filename.text)))
            dataMeshPaths.append(structDataMeshPaths)
    else:               # NEW AND CORRECT XML CONVENTION
        subject = dataSetTree.find('subject')
        visits = subject.findall('visit')
        for visit in visits:
            dataAges.append(float(visit.find('age').text))
        for struct in structureNames:
            structDataMeshPaths = []
            for visit in visits:
                filenames = visit.findall('filename')
                for filename in filenames:
                    objectId = filename.attrib['object_id']
                    if objectId == struct:
                        structDataMeshPaths.append(os.path.normpath(os.path.join(pathToRegressionFolder, "output", filename.text)))
            dataMeshPaths.append(structDataMeshPaths)
    return dataAges, dataMeshPaths

def computeVolume(meshPath):
    #Reading
    reader = vtkPolyDataReader()
    reader.SetFileName(meshPath)
    reader.Update()
    #First volume
    mass = vtkMassProperties()
    mass.SetInputConnection(reader.GetOutputPort())
    mass.Update()
    volume = mass.GetVolume()
    assert 0 <= volume, "Warning : negative volume, that's weird."
    return volume

def computeRegressionVolumes(pathToRegressionFolder, structureNames, t0, tn, nbTimePoints):
    vectorRegressionAges = []
    arrayRegressionVolumes = []
    for tp in range(nbTimePoints):
        vectorRegressionAges.append(t0 + tp * (tn - t0)/(nbTimePoints - 1))
    outputFiles = os.listdir(os.path.join(pathToRegressionFolder, "output"))
    for struct in structureNames:
        structRegressionVolumes = [None] * nbTimePoints
        for filename in outputFiles:
            if filename.find(struct + "_trajectory___t_")>=0:
                index = extractIndex(filename)
                vol = computeVolume(os.path.join(pathToRegressionFolder, "output", filename))
                structRegressionVolumes[index] = vol
        arrayRegressionVolumes.append(structRegressionVolumes)
    return vectorRegressionAges, arrayRegressionVolumes

def computeObservedVolumes(pathToRegressionFolder, structureNames, dataAges, dataMeshPaths):
    vectorObservationAges = dataAges
    arrayObservedVolumes = []
    structIndex = 0
    for struct in structureNames:
        structObservedVolumes = []
        for filename in dataMeshPaths[structIndex]:
            vol = computeVolume(os.path.join(pathToRegressionFolder, "output", filename))
            structObservedVolumes.append(vol)
        arrayObservedVolumes.append(structObservedVolumes)
        structIndex = structIndex + 1
    return vectorObservationAges, arrayObservedVolumes

def computeTotalVolumes(structureNames, vectorRegressionAges, arrayRegressionVolumes, vectorObservationAges, arrayObservedVolumes):
    vectorTotalRegressionVolumes = []
    for ageIndex in range(len(vectorRegressionAges)):
        avg = 0.0
        for structIndex in range(len(structureNames)):
            avg = avg + arrayRegressionVolumes[structIndex][ageIndex]
        vectorTotalRegressionVolumes.append(avg)
    vectorTotalObservedVolumes = []
    for ageIndex in range(len(vectorObservationAges)):
        avg = 0.0
        for structIndex in range(len(structureNames)):
            avg = avg + arrayObservedVolumes[structIndex][ageIndex]
        vectorTotalObservedVolumes.append(avg)
    return vectorTotalRegressionVolumes, vectorTotalObservedVolumes

def plotVolumes(pathToRegressionFolder, structureNames, vectorRegressionAges, arrayRegressionVolumes, vectorTotalRegressionVolumes, vectorObservationAges, arrayObservedVolumes, vectorTotalObservedVolumes):
    plt.clf()
    plt.plot(vectorObservationAges, vectorTotalObservedVolumes, 'k', label="raw data total volume evolution")
    plt.plot(vectorRegressionAges, vectorTotalRegressionVolumes, 'k--', label="reg data total volume evolution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pathToRegressionFolder, "output", "volumeEvolution_total.pdf"))
    # plt.show()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    plt.clf()
    for structIndex in range(len(structureNames)):
        plt.plot(vectorObservationAges, arrayObservedVolumes[structIndex], colors[structIndex], label="raw (solid) and reg (dashed) " + structureNames[structIndex])
        plt.plot(vectorRegressionAges, arrayRegressionVolumes[structIndex], colors[structIndex]+'--')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pathToRegressionFolder, "output", "volumeEvolution_detailed.pdf"))
    plt.clf()
    # plt.show()

assert len(sys.argv)>=2, "Usage: plotVolumeEvolutionForRegression.py pathToRegressionFolder"
pathToRegressionFolder = sys.argv[1]

# Main.
structureNames, t0, tn, nbTimePoints = extractInfoFromModel(pathToRegressionFolder)
dataAges, dataMeshPaths = extractInfoFromDataSet(pathToRegressionFolder, structureNames)
vectorRegressionAges, arrayRegressionVolumes = computeRegressionVolumes(pathToRegressionFolder, structureNames, t0, tn, nbTimePoints)
vectorObservationAges, arrayObservedVolumes = computeObservedVolumes(pathToRegressionFolder, structureNames, dataAges, dataMeshPaths)
vectorTotalRegressionVolumes, vectorTotalObservedVolumes = computeTotalVolumes(structureNames, vectorRegressionAges, arrayRegressionVolumes, vectorObservationAges, arrayObservedVolumes)
plotVolumes(pathToRegressionFolder, structureNames, vectorRegressionAges, arrayRegressionVolumes, vectorTotalRegressionVolumes, vectorObservationAges, arrayObservedVolumes, vectorTotalObservedVolumes)
