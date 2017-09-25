import random as rd
import vtk
from vtk import vtkPolyDataReader, vtkBooleanOperationPolyDataFilter, vtkMassProperties, vtkPolyDataWriter, vtkIntersectionPolyDataFilter
import numpy as np
import sys


##Warning : the dice distance function below is somewhat unstable when meshes are highly non-convex. It gives accurate results for Hippocampus, but expect weird behaviors for the cortex...

def computeDiceDistance(mesh1Path, mesh2Path):
    if mesh1Path==mesh2Path:
        return 1.
    #Reading
    reader1 = vtkPolyDataReader()
    reader1.SetFileName(mesh1Path)
    reader1.Update()
    reader2 = vtkPolyDataReader()
    reader2.SetFileName(mesh2Path)
    reader2.Update()
    #First volume
    mass1 = vtkMassProperties()
    mass1.SetInputConnection(reader1.GetOutputPort())
    mass1.Update()
    volume1 = mass1.GetVolume()
    #Second volume
    mass2 = vtkMassProperties()
    mass2.SetInputConnection(reader2.GetOutputPort())
    mass2.Update()
    volume2 = mass2.GetVolume()
    #Intersection
    # intersectionOperation = vtkIntersectionPolyDataFilter()
    intersectionOperation = vtkBooleanOperationPolyDataFilter()
    intersectionOperation.SetOperationToIntersection()
    intersectionOperation.SetInputConnection(0, reader1.GetOutputPort())
    intersectionOperation.SetInputConnection(1, reader2.GetOutputPort())
    intersectionOperation.Update()
    #Volume of the intersection
    massInter = vtkMassProperties()
    massInter.SetInputConnection(intersectionOperation.GetOutputPort())
    massInter.Update()
    intersectionVolume = massInter.GetVolume()
    dice = 2*intersectionVolume/(volume1+volume2)
    assert 0 <= dice <= 1, "Warning : suspicious behavior detected in the intersection filter."
    return dice


assert len(sys.argv)>=3, "Usage: reorientMesh.py mesh1Path mesh2Path"
meshPath1 = sys.argv[1]
meshPath2 = sys.argv[2]


print(computeDiceDistance(meshPath1, meshPath2))
