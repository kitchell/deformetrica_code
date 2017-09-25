import random as rd
import vtk
from vtk import vtkPolyDataWriter, vtkPolyDataReader,  vtkTransform, vtkMatrix4x4, vtkTransformPolyDataFilter
import numpy as np
import sys

assert len(sys.argv)>3, "Usage: reorientMesh.py meshPath transformationMatrixPath savePath"

meshPath = sys.argv[1]
transformationMatrixPath = sys.argv[2]
saveName = sys.argv[3]

#Setup a reader for the mesh (assuming it's polydata)
reader = vtk.vtkPolyDataReader()
reader.SetFileName(meshPath)
reader.Update()

#Read the matrix file:
transformationMatrix = np.loadtxt(transformationMatrixPath)
x,y = transformationMatrix.shape
print(x,y)
assert x==3 and y==4, "Wrong dimension of input matrix : %i %i, expecting 4x4" % x % y

transformationMatrix.transpose()

#Copy this array into a vtkMatrix.
vtkTransformationMatrix = vtkMatrix4x4()
for i in range(0,3):
    for j in range(0,4):
        vtkTransformationMatrix.SetElement(i,j,transformationMatrix[i,j])

lastLine = np.array([0.,0.,0.,1.])
for j in range(4):
    vtkTransformationMatrix.SetElement(3,j,lastLine[j])


#Creates a linear transformation from this matrix object.
linearTransform = vtkTransform()
linearTransform.SetMatrix(vtkTransformationMatrix)
#print(linearTransform)

#Define a filter to apply the transformation to the output of the reader.
transformFilter = vtkTransformPolyDataFilter()
transformFilter.SetTransform(linearTransform)
transformFilter.SetInputConnection(reader.GetOutputPort())
transformFilter.Update()

#Write the output of the filter into a new vtk file.
writer = vtkPolyDataWriter()
writer.SetInputConnection(transformFilter.GetOutputPort())
writer.SetFileName(saveName)
writer.Update()
