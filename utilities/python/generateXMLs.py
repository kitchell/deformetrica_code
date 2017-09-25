import numpy as np
import sys
import os
import xml.etree.ElementTree as ET

#the xmls generated with these functions are unindented, but most text editors know how to indent xml files.

def generateModelXml(t0, tn, numberOfTimePoints, deformationKernelWidth, deformationKernelType, templateMeshes, templateIds, objectTypes, structureKernelTypes, savePath):
    """Generate a model xml file from input parameters :
    -numberOfTimePoints
    -deformationKernelWidth and deformationKernelType : parameters of the deformation
    -templateMeshes, templateIds, objectTypes, structureKernelTypes : descriptions of the template objects
    -savePath : path used to save the generated xml file.
    """
    a = ET.Element('model')
    modelType = ET.SubElement(a,'model-type')
    modelType.text = "Regression"
    #Deformation-parameters
    deformationParameters = ET.SubElement(a,'deformation-parameters')
    t0 = ET.SubElement(deformationParameters,'t0')
    t0.text = t0
    tn = ET.SubElement(deformationParameters,'tn')
    tn.text = tn
    kernelWidth = ET.SubElement(deformationParameters,'kernel-width')
    kernelWidth.text = deformationKernelWidth
    kernelType = ET.SubElement(deformationParameters,'kernel-type')
    kernelType.text = deformationKernelType
    numberOfTimePoints = ET.SubElement(deformationParameters,'number-of-timepoints')
    numberOfTimePoints.text = numberOfTimePoints
    #Template
    template = ET.SubElement(a,'template')
    for i,path in enumerate(templateMeshes):
        obj = ET.SubElement(template, 'object')
        obj.set('id', templateIds[i])
        objectType = ET.SubElement(obj, 'deformable-object-type')
        objectType.text = objectTypes[i]
        objKernelType = ET.SubElement(obj, 'kernel-type')
        objKernelType.text = structureKernelTypes[i]
        objName = ET.SubElement(obj, 'filename')
        objName.text = path
    doc = ET.ElementTree(a)
    doc.write(savePath)


def generateDataSetXML(meshes, meshesId, savePath):
    """Generate a data_set xml file from input parameters :
   -meshesId : ids of the meshes,
   -meshesTime : ages (useful in a regression)
   -meshesPath : paths of the meshes or images
   """
    a = ET.Element('data-set')
    sub2 = ET.SubElement(a, 'subject')
    sub2.set('id', "singleSubjectHere")
    visit2 = ET.SubElement(sub2,'visit')
    for i,mesh in enumerate(meshes):
        f = ET.SubElement(visit2, 'filename')
        f.text = mesh
        f.set('object_id', meshesId[i])
    doc = ET.ElementTree(a)
    doc.write(savePath)


def generateOptimizationParameterXml(freezeTemplate, savePath):
    a = ET.Element('optimization-parameters')
    freezeT = ET.Subelement(a, 'freeze-template')
    freezeT.text = freezeTemplate
    doc = ET.ElementTree(a)
    doc.write(savePath)

