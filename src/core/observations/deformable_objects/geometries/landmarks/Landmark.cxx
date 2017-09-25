/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "Landmark.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTransformFilter.h"

#include <cstring>
#include <iostream>
#include <sstream>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
Landmark<ScalarType, Dimension>
::Landmark() : Superclass() {
  this->SetLandmarkType();
  m_NumberOfPoints = 0;
}

template<class ScalarType, unsigned int Dimension>
Landmark<ScalarType, Dimension>
::~Landmark() {

}

template<class ScalarType, unsigned int Dimension>
Landmark<ScalarType, Dimension>
::Landmark(const Landmark &other) : Superclass(other) {
  this->SetLandmarkType();

  if (other.IsOfUnstructuredKind()) {
    ///Get a vtkUnstructuredGrid pointer from the vtkPointSet one.
    vtkSmartPointer<vtkUnstructuredGrid> realPointer = vtkUnstructuredGrid::SafeDownCast(other.m_PointSet);

    ///Deep copy the data into a new pointer.
    vtkSmartPointer<vtkUnstructuredGrid> copy = vtkSmartPointer<vtkUnstructuredGrid>::New();
    copy->DeepCopy(realPointer);

    ///Upcast the copy and set the m_PointSet attribute.
    vtkSmartPointer<vtkPointSet> data = vtkSmartPointer<vtkPointSet>(copy);
    m_PointSet = data;
  } else {
    ///Get a vtkPolyData pointer from the vtkPointSet one.
    vtkSmartPointer<vtkPolyData> realPointer = vtkPolyData::SafeDownCast(other.m_PointSet);

    ///Deepcopy the vtkPolyData into a new pointer.
    vtkSmartPointer<vtkPolyData> copy = vtkSmartPointer<vtkPolyData>::New();
    copy->DeepCopy(realPointer);

    //Upcast the copy and set the m_PointSet attribute.
    vtkSmartPointer<vtkPointSet> data = vtkSmartPointer<vtkPointSet>(copy);
    m_PointSet = data;
  }

  m_PointCoordinates = other.m_PointCoordinates;

  m_NumberOfPoints = other.m_NumberOfPoints;

}

template<class ScalarType, unsigned int Dimension>
Landmark<ScalarType, Dimension>
::Landmark(const Landmark &example, const MatrixType &LandmarkPoints) : Superclass(example, LandmarkPoints) {
  this->SetLandmarkType();

  m_NumberOfPoints = example.m_NumberOfPoints;

  if (m_NumberOfPoints != LandmarkPoints.rows())
    throw std::runtime_error("Number of LandmarkPoints mismatch in copy/update Landmark constructor");

  if (example.IsOfUnstructuredKind()) {
    vtkSmartPointer<vtkUnstructuredGrid> copy = vtkSmartPointer<vtkUnstructuredGrid>::New();
    copy->DeepCopy(example.m_PointSet);
    m_VTKMutex.Lock();
    m_PointSet = vtkSmartPointer<vtkPointSet>(copy);
    m_VTKMutex.Unlock();
  } else {
    vtkSmartPointer<vtkPolyData> copy = vtkSmartPointer<vtkPolyData>::New();
    copy->DeepCopy(example.m_PointSet);
    m_VTKMutex.Lock();
    m_PointSet = vtkSmartPointer<vtkPointSet>(copy);
    m_VTKMutex.Unlock();
  }

  this->UpdatePolyDataWithPointCoordinates(LandmarkPoints);

  this->Update();
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Landmark<ScalarType, Dimension>
::SetPolyData(vtkPointSet *pointSet) {
  if (Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel().compare("LPS") != 0) {
    std::cout << "Reorienting polydata from "
              << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel() << " to LPS" << std::endl;
    m_PointSet = ReorientPolyData(pointSet, false);
  } else
    m_PointSet = vtkSmartPointer<vtkPointSet>(pointSet);

  m_NumberOfPoints = m_PointSet->GetNumberOfPoints();

  m_PointCoordinates.set_size(m_NumberOfPoints, Dimension);

  for (unsigned int i = 0; i < m_NumberOfPoints; i++) {
    // Here we used 3 since 2D points still have a z-coordinate that is equal to 0 in vtkPolyData.
    // This coordinate is removed in m_WorkingPointCoordinates
    double p[3];
    m_VTKMutex.Lock();
    m_PointSet->GetPoint(i, p);
    m_VTKMutex.Unlock();
    for (int dim = 0; dim < Dimension; dim++)
      m_PointCoordinates(i, dim) = p[dim];
  }

  this->SetModified();
}

template<class ScalarType, unsigned int Dimension>
vtkSmartPointer<vtkPointSet>
Landmark<ScalarType, Dimension>
::ReorientPolyData(vtkPointSet *pointSet, bool inverse) const {
  double trfmtx[16];
  MatrixType orient_mtx;

  if (inverse == false)
    orient_mtx = Superclass::m_AnatomicalOrientation.GetChangeOfBasisMatrix();
  else
    orient_mtx = Superclass::m_AnatomicalOrientation.GetInverseChangeOfBasisMatrix();

  int idx = 0;
  for (unsigned int i = 0; i < Dimension; i++) {
    for (unsigned int j = 0; j < Dimension; j++) {
      trfmtx[idx++] = orient_mtx(i, j);
      std::cout << trfmtx[idx - 1] << " ";
    }
    trfmtx[idx++] = 0;
    std::cout << trfmtx[idx - 1] << std::endl;
  }
  for (unsigned int j = 0; j < Dimension; j++) {
    trfmtx[idx++] = 0;
    std::cout << trfmtx[idx - 1] << " ";
  }
  trfmtx[idx++] = 1;
  std::cout << trfmtx[idx - 1] << std::endl;

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix(trfmtx);

  ///The way we reorient the data depends on the actual type of data : vtkPolyData or vtkUnstructuredGrid.
  if (this->IsOfUnstructuredKind()) {
    //Here the data is a vtkUnstructuredGrid.
    vtkSmartPointer<vtkTransformFilter> transformUnstructuredData = vtkSmartPointer<vtkTransformFilter>::New();

    //Downcast the vtkPointSet to a vtkUnstructuredGrid.
    vtkSmartPointer<vtkUnstructuredGrid> downcastedPointer = vtkUnstructuredGrid::SafeDownCast(pointSet);

#if VTK_MAJOR_VERSION <= 5
    transformUnstructuredData->SetInput(toModify);
    //transformUnstructuredData->SetInputConnection(polyData->GetProducerPort());
#else
    transformUnstructuredData->SetInputData(downcastedPointer);
#endif
    transformUnstructuredData->SetTransform(transform);
    transformUnstructuredData->Update();
    vtkSmartPointer<vtkPointSet> pointSetToReturn = transformUnstructuredData->GetOutput();
    pointSetToReturn->SetPoints(pointSetToReturn->GetPoints());//This step is necessary because the transform only returns a vtkPointSet and not an unstructuredgrid.
    return pointSetToReturn;

  } else {
    //Here the data is a vtkPolyData.
    vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyData = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

    //Downcast the vtkPointSet to a vtkPolyData.
    vtkSmartPointer<vtkPolyData> downcastedPointer = vtkPolyData::SafeDownCast(pointSet);

    vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();

#if VTK_MAJOR_VERSION <= 5
    transformPolyData->SetInput(polyData);
    //transformPolyData->SetInputConnection(polyData->GetProducerPort());
#else
    transformPolyData->SetInputData(downcastedPointer);
#endif
    transformPolyData->SetTransform(transform);
    transformPolyData->SetOutput(newPolyData);
    transformPolyData->Update();

    return transformPolyData->GetOutput();
  }
}

template<class ScalarType, unsigned int Dimension>
void
Landmark<ScalarType, Dimension>
::UpdatePolyDataWithPointCoordinates(const MatrixType &Y) {
  // A polydata needs to be set before updating Point Coordinates
  if (m_PointSet == NULL)
    throw std::runtime_error("a VTK PolyData should be set before setting new point coordinates");

  if (Y.rows() != m_PointSet->GetNumberOfPoints())
    throw std::runtime_error("number of points mismatched");
  if (Y.columns() != Dimension)
    throw std::runtime_error("Dimension mismatched");

  for (unsigned int i = 0; i < m_NumberOfPoints; i++) {
    double p[3];
    p[2] = 0.0;
    for (int dim = 0; dim < Dimension; dim++)
      p[dim] = Y(i, dim);

    m_VTKMutex.Lock();
    m_PointSet->GetPoints()->SetPoint(i, p);
    m_VTKMutex.Unlock();
  }

  m_PointCoordinates = Y;

  this->SetModified();
}

template<class ScalarType, unsigned int Dimension>
int Landmark<ScalarType, Dimension>::GetCellForPosition(const VectorType &pos, double *weights) {
  ///Declaring the parameters for FindCell.
  int subId = -1;
  double x[3];
  x[2] = 0.;
  for (int i = 0; i < Dimension; i++)
    x[i] = pos[i];
  double pcoords[3];
  this->m_VTKMutex.Lock();
  vtkIdType cellId = m_PointSet->FindCell(x, NULL, -1, 1e-1, subId, pcoords, weights);//TODO withdraw this 1e-5
  this->m_VTKMutex.Unlock();

  return cellId;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Landmark<ScalarType, Dimension>
::Update() {
  if (m_PointSet == NULL)
    throw std::runtime_error("A VTK PolyData should have been set in Landmark class or its children classes");

  if (this->IsModified()) {
    this->UpdateBoundingBox();
  }

  this->UnSetModified();
}

template<class ScalarType, unsigned int Dimension>
ScalarType
Landmark<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<Superclass> target) {
  if (this->GetType() != target->GetType())
    throw std::runtime_error("Abstract geometries types mismatched");

  const std::shared_ptr<const Landmark> targetLandmark = std::static_pointer_cast<const Landmark>(target);

  if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
    throw std::runtime_error("Landmark object should have the same number of points");

  MatrixType D = this->GetPointCoordinates() - targetLandmark->GetPointCoordinates();

  ScalarType match = D.frobenius_norm();
  match *= match;

  return match;
}

template<class ScalarType, unsigned int Dimension>
MatrixType
Landmark<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<Superclass> target) {
  if (this->GetType() != target->GetType())
    throw std::runtime_error("Deformable objects types mismatched");

  const std::shared_ptr<const Landmark> targetLandmark = std::static_pointer_cast<const Landmark>(target);

  if (m_NumberOfPoints != targetLandmark->GetPointCoordinates().rows())
    throw std::runtime_error("Landmarks object should have the same number of points");

  MatrixType gradMatch = this->GetPointCoordinates() - targetLandmark->GetPointCoordinates();
  gradMatch *= 2.0f;

  return gradMatch;
}

template<class ScalarType, unsigned int Dimension>
void
Landmark<ScalarType, Dimension>
::WriteObject(std::string str) const {
  vtkSmartPointer<vtkPointSet> outData = m_PointSet;

  if (Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel().compare("LPS") != 0) {
    std::cout << "Reorienting polydata to "
              << Superclass::m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel()
              << " before writing output file" << endl;

    outData = ReorientPolyData(m_PointSet, true);
  }

  if (this->IsOfUnstructuredKind()) {
    std::cout << "Writing object as a vtkUnstructuredGrid " << std::endl;
    vtkSmartPointer<vtkUnstructuredGrid> downcastedPointer = vtkUnstructuredGrid::SafeDownCast(outData);
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(str.c_str());
#if (VTK_MAJOR_VERSION == 5)
    writer->SetInput(downcastedPointer);
#else
    writer->SetInputData(downcastedPointer); //
#endif
    writer->Update();

  } else {
    vtkSmartPointer<vtkPolyData> downcastedPointer = vtkPolyData::SafeDownCast(outData);
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(str.c_str());
#if (VTK_MAJOR_VERSION == 5)
    writer->SetInput(downcastedPointer);
#else
    writer->SetInputData(downcastedPointer); //
#endif
    writer->Update();
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Landmark<ScalarType, Dimension>
::UpdateBoundingBox() {

  VectorType Min(Dimension, 1e20);
  VectorType Max(Dimension, -1e20);

  for (int i = 0; i < m_NumberOfPoints; i++) {
    VectorType p = m_PointCoordinates.get_row(i);
    for (unsigned int dim = 0; dim < Dimension; dim++) {
      Min[dim] = (p[dim] < Min[dim]) ? p[dim] : Min[dim];
      Max[dim] = (p[dim] > Max[dim]) ? p[dim] : Max[dim];
    }
  }

  Superclass::m_BoundingBox.set_column(0, Min);
  Superclass::m_BoundingBox.set_column(1, Max);
}

template class Landmark<double,2>;
template class Landmark<double,3>;
