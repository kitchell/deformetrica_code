
import random as rd
import vtk
from vtk import vtkPolyDataWriter, vtkDelaunay3D, vtkUnstructuredGridWriter, vtkMath, vtkDelaunay2D.
import numpy as np


#Given three lengths, generates a set of uniformely spread points within the ellipsoid.
def generateEllipsoidPoints(a1, a2, a3, numPoints):
    l = []
    i = 0
    while (i < numPoints):
        x = rd.uniform(-a1, a1)
        y = rd.uniform(-a2, a2)
        z = rd.uniform(-a3, a3)
        if ((x/a1)**2 + (y/a2)**2 + (z/a3)**2) < 1:#in the ellipsoid
            l.append([i,x,y,z])
            i += 1
    return l

#Generates a set of points uniformely scattered within a cube of size 1.
def generateCubePoints(numPoints):
    l = []
    l.append([0,1,1,1])
    l.append([1,1,1,-1])
    l.append([2,1,-1,1])
    l.append([3,1,-1,-1])
    l.append([4,-1,1,1])
    l.append([5,-1,1,-1])
    l.append([6,-1,-1,1])
    l.append([7,-1,-1,-1])
    i = 8
    while (i<numPoints):
        x = rd.uniform(-1,1)
        y = rd.uniform(-1,1)
        z = rd.uniform(-1,1)
        l.append([i,x,y,z])
        i+=1
    return l

#Given a set of points, save the corresponding full volume mesh into an unstructured grid file.
def generateVTKFromPoints(listOfPoints):
    math = vtk.vtkMath()
    points = vtk.vtkPoints()
    for elt in listOfPoints:
        (i,x,y,z) = elt
        points.InsertPoint(i,x,y,z)

    profile = vtk.vtkPolyData()
    profile.SetPoints(points)

    delny = vtk.vtkDelaunay3D()
    delny.SetInputData(profile)
    delny.SetTolerance(0.0001)
    delny.SetAlpha(0)
    print("alpha :",delny.GetAlpha())

    print("Output :", delny.GetOutput())
    delny.BoundingTriangulationOff()

    toSave = vtk.vtkUnstructuredGridWriter()
    toSave.SetInputConnection(delny.GetOutputPort())

    toSave.SetFileName("Cube"+str(len(listOfPoints))+".vtk")
    toSave.Write()

#Given a list of 2D points, generates de corresponding vtkPolyData object abd saves it.
def generateVTK2DFromPoints(listOfPoints, name):
    math = vtk.vtkMath()
    points = vtk.vtkPoints()
    for elt in listOfPoints:
        (i,x,y) = elt
        points.InsertPoint(i,x,y,0)

    profile = vtk.vtkPolyData()
    profile.SetPoints(points)

    delny = vtk.vtkDelaunay2D()
    delny.SetInputData(profile)
    delny.SetTolerance(0.0001)
    delny.SetAlpha(0)
    print("alpha :",delny.GetAlpha())
    print("Output :", delny.GetOutput())
    delny.BoundingTriangulationOff()

    toSave = vtk.vtkPolyDataWriter()
    toSave.SetInputConnection(delny.GetOutputPort())

    toSave.SetFileName(name)
    toSave.Write()

#Generate a set of 2D points scattered over a square.
def generateSquarePoints(nbPoints):
        l = []
        l.append([0,1,1])
        l.append([1,1,-1])
        l.append([2,1,-1])
        l.append([3,-1,-1])
        i = 4
        while (i<nbPoints):
            x = rd.uniform(-1,1)
            y = rd.uniform(-1,1)
            l.append([i,x,y])
            i+=1
        return l

#Geenrate a set of 2D points scattered over the unit disk.
def generateDiskPoints(nbPoints):
        l = []
        i=0
        while (i<nbPoints):
            x = rd.uniform(-1,1)
            y = rd.uniform(-1,1)
            if (x**2+y**2<1):
                l.append([i,x,y])
                i+=1
        return l

