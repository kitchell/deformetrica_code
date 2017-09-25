/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "OrientedVolumeMesh.h"

#include "KernelFactory.h"


#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkTriangleFilter.h"
#include "myvtkPolyDataNormals.h"

#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
OrientedVolumeMesh<ScalarType, Dimension>
::OrientedVolumeMesh() : Superclass()
{
  this->SetOrientedVolumeMeshType();
  m_KernelWidth = 0;
  m_KernelType = null;
  m_Reorient = false;
}



template <class ScalarType, unsigned int Dimension>
OrientedVolumeMesh<ScalarType, Dimension>
::~OrientedVolumeMesh()
{}



template <class ScalarType, unsigned int Dimension>
OrientedVolumeMesh<ScalarType, Dimension>
::OrientedVolumeMesh(const OrientedVolumeMesh& other) : Superclass(other)
{
  this->SetOrientedVolumeMeshType();

  m_Centers = other.m_Centers;
  m_Volumes = other.m_Volumes;
  m_NumCells = other.m_NumCells;

  m_KernelWidth = other.m_KernelWidth;
  m_KernelType = other.m_KernelType;
  m_NormSquared = other.m_NormSquared;

  m_Reorient = other.m_Reorient;
}



template <class ScalarType, unsigned int Dimension>
OrientedVolumeMesh<ScalarType, Dimension>
::OrientedVolumeMesh(const OrientedVolumeMesh& ex, const MatrixType& LP) : Superclass(ex, LP)
{
  this->SetOrientedVolumeMeshType();

  m_KernelWidth = ex.m_KernelWidth;
  m_KernelType = ex.m_KernelType;

  m_Reorient = ex.m_Reorient;

  // this is required since the call to Superclass(ex, LP) sets m_IsModified to false
  this->SetModified();
  this->Update();
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
OrientedVolumeMesh<ScalarType, Dimension>
::Update()
{
  if (this->IsModified())
  {
    Superclass::Update();
    this->CheckMeshAndNormals();
    this->UpdateCentersNormals();
    // replace the BoundingBox computed with vertices by the bounding box computed with centers of triangles
    this->UpdateBoundingBox();
    this->UpdateSelfNorm();
  }
  this->UnSetModified();
}



template<class ScalarType, unsigned int Dimension>
ScalarType OrientedVolumeMesh<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<AbstractGeometryType> target)
{

  if (this->GetType() != target->GetType())
    throw std::runtime_error("Abstract Geometries types mismatched");

  const std::shared_ptr<const OrientedVolumeMesh> targetOrientedVolumeMesh
      = std::static_pointer_cast<const OrientedVolumeMesh>(target);

  if (m_KernelWidth != targetOrientedVolumeMesh->GetKernelWidth())
    throw std::runtime_error("Kernel width of surface currents mismatched");

  this->Update();
  // target->Update();

  MatrixType targCenters = targetOrientedVolumeMesh->GetCenters();
  MatrixType targVolumes = targetOrientedVolumeMesh->GetVolumes();

  ScalarType match = targetOrientedVolumeMesh->GetNormSquared() + this->GetNormSquared();

  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

  MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

  std::shared_ptr<KernelType> kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
  kernelObject->SetKernelWidth(m_KernelWidth);
  kernelObject->SetSources(m_Centers);
  kernelObject->SetWeights(m_Volumes);

  MatrixType SdotT = kernelObject->Convolve(targCenters);
  for (int i = 0; i < targetOrientedVolumeMesh->GetNumberOfCells(); i++)
    match -= 2.0f * dot_product(SdotT.get_row(i), targVolumes.get_row(i));

  return match;
}



template<class ScalarType, unsigned int Dimension>
MatrixType OrientedVolumeMesh<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target)
{
  if (this->GetType() != target->GetType())
    throw std::runtime_error("Abstract Geometries types mismatched");

  const std::shared_ptr<const OrientedVolumeMesh> targetOrientedVolumeMesh =
      std::static_pointer_cast<const OrientedVolumeMesh>(target);

  if (m_KernelWidth != targetOrientedVolumeMesh->GetKernelWidth())
    throw std::runtime_error("Kernel width of surface currents mismatched");

  this->Update();
  // target->Update();

  MatrixType targCenters = targetOrientedVolumeMesh->GetCenters();
  MatrixType targVolumes = targetOrientedVolumeMesh->GetVolumes();

  MatrixType Pts = this->GetPointCoordinates();

  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

  MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

  std::shared_ptr<KernelType> kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
  kernelObject->SetKernelWidth(m_KernelWidth);

  kernelObject->SetSources(m_Centers);
  kernelObject->SetWeights(m_Volumes);
  VectorType KtauS(kernelObject->Convolve(m_Centers).get_column(0));
  MatrixType gradKtauS = kernelObject->ConvolveGradient(m_Centers, m_Volumes);

  kernelObject->SetSources(targCenters);
  kernelObject->SetWeights(targVolumes);
  VectorType KtauT(kernelObject->Convolve(m_Centers).get_column(0));
  MatrixType gradKtauT = kernelObject->ConvolveGradient(m_Centers, m_Volumes);

  //This is a third or a fourth coming from the contribution of a vertex to the center divided by 2
  float dimensionFactor = (Dimension+1.0f)/2.0f;

  MatrixType gradKtau = (gradKtauS - gradKtauT) / dimensionFactor;
  MatrixType gradmatch(this->GetNumberOfPoints(), Dimension, 0.0);

  //This if condition is disgusting, but not computationally expensive :) (could be replaced with a general chain rule computation)
  if (Dimension==3)
  {
    for (int f = 0; f < m_NumCells; f++) {
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

      m_VTKMutex.Lock();
      Superclass::m_PointSet->GetCellPoints(f, ptIds);
      m_VTKMutex.Unlock();

      int ind0 = ptIds->GetId(0);
      int ind1 = ptIds->GetId(1);
      int ind2 = ptIds->GetId(2);
      int ind3 = ptIds->GetId(3);

      //The variation of the normal caused by each of the point of the tetrahedron.
      VectorType v0 = cross_3d(Pts.get_row(ind3) - Pts.get_row(ind1), Pts.get_row(ind2) - Pts.get_row(ind1));
      VectorType v1 = cross_3d(Pts.get_row(ind2) - Pts.get_row(ind0), Pts.get_row(ind3) - Pts.get_row(ind0));
      VectorType v2 = cross_3d(Pts.get_row(ind3) - Pts.get_row(ind0), Pts.get_row(ind1) - Pts.get_row(ind0));
      VectorType v3 = cross_3d(Pts.get_row(ind1) - Pts.get_row(ind0), Pts.get_row(ind2) - Pts.get_row(ind0));

      //The variation of the kernel when a point moves
      ScalarType Ktau = KtauS(f) - KtauT(f);

      gradmatch.set_row(ind0, gradmatch.get_row(ind0) + 1.0f/3.0f * v0 * Ktau);
      gradmatch.set_row(ind1, gradmatch.get_row(ind1) + 1.0f/3.0f * v1 * Ktau);
      gradmatch.set_row(ind2, gradmatch.get_row(ind2) + 1.0f/3.0f * v2 * Ktau);
      gradmatch.set_row(ind3, gradmatch.get_row(ind3) + 1.0f/3.0f * v3 * Ktau);

      gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau.get_row(f));
      gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau.get_row(f));
      gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau.get_row(f));
      gradmatch.set_row(ind3, gradmatch.get_row(ind3) + gradKtau.get_row(f));
    }
  }

  else
  {
    for (int f = 0; f < m_NumCells; f++) {
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

      m_VTKMutex.Lock();
      Superclass::m_PointSet->GetCellPoints(f, ptIds);
      m_VTKMutex.Unlock();

      int ind0 = ptIds->GetId(0);
      int ind1 = ptIds->GetId(1);
      int ind2 = ptIds->GetId(2);

      VectorType point0 = Pts.get_row(ind0);
      VectorType point1 = Pts.get_row(ind1);
      VectorType point2 = Pts.get_row(ind2);

      //The variation of the volume caused by each of the point of the tetrahedron.
      VectorType v0(2), v1(2), v2(2);
      v0(0) = point1(1) - point2(1);
      v0(1) = point2(0) - point1(0);
      v1(0) = point2(1) - point0(1);
      v1(1) = point0(0) - point2(0);
      v2(0) = point0(1) - point1(1);
      v2(1) = point1(0) - point0(0);

      //The variation of the kernel when a point moves
      ScalarType Ktau = KtauS(f) - KtauT(f);

      gradmatch.set_row(ind0, gradmatch.get_row(ind0) + v0 * Ktau);
      gradmatch.set_row(ind1, gradmatch.get_row(ind1) + v1 * Ktau);
      gradmatch.set_row(ind2, gradmatch.get_row(ind2) + v2 * Ktau);

      gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau.get_row(f));
      gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau.get_row(f));
      gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau.get_row(f));
    }
  }

  return gradmatch;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
OrientedVolumeMesh<ScalarType, Dimension>
::CheckMeshAndNormals()
{
  m_NumCells = Superclass::m_PointSet -> GetNumberOfCells();
}



template <class ScalarType, unsigned int Dimension>
void
OrientedVolumeMesh<ScalarType, Dimension>
::UpdateCentersNormals()
{
  const MatrixType Pts = this->GetPointCoordinates();
  m_Centers.set_size(m_NumCells, Dimension);
  m_Volumes.set_size(m_NumCells, 1);

  //The number of triangles in a parallelogram or of tetrahedron in a cube (to divide the determinant).
  float dimensionFactor = (Dimension==3) ? (6.0f) : (2.0f);

  for (unsigned int i = 0; i < m_NumCells; i++)
  {
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    m_VTKMutex.Lock();
    Superclass::m_PointSet->GetCellPoints(i, ptIds);
    m_VTKMutex.Unlock();

    if (ptIds->GetNumberOfIds() != 4 && ptIds->GetNumberOfIds() != 3)
      throw std::runtime_error("Not an admissible cell!");

    MatrixType edges(Dimension, Dimension);
    VectorType center = Pts.get_row(ptIds->GetId(0));

    for (int d=1; d<Dimension+1; d++)
    {
      center += Pts.get_row(ptIds->GetId(d));
      edges.set_row(d-1, Pts.get_row(ptIds->GetId(d)) - Pts.get_row(ptIds->GetId(0)));
    }
    center /= (Dimension+1.0f);

    m_Centers.set_row(i, center);

    ScalarType x = std::abs(edges.determinant() / dimensionFactor); //this is the wedge product of three edges from the same vertice.

    m_Volumes.set_row(i, x);

  }
}



template <class ScalarType, unsigned int Dimension>
void
OrientedVolumeMesh<ScalarType, Dimension>
::UpdateSelfNorm()
{
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

  std::shared_ptr<KernelType> kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), this->GetBoundingBox());

  kernelObject->SetKernelWidth(m_KernelWidth);
  kernelObject->SetSources(m_Centers);
  kernelObject->SetWeights(m_Volumes);

  MatrixType selfKW = kernelObject->Convolve(m_Centers);

  m_NormSquared = 0;

  for (unsigned int i = 0; i < m_NumCells; i++)
    m_NormSquared += dot_product(selfKW.get_row(i), m_Volumes.get_row(i));
}



template <class ScalarType, unsigned int Dimension>
void
OrientedVolumeMesh<ScalarType, Dimension>
::UpdateBoundingBox()
{
  // overload bounding box using the centers instead of the points (useful when using the results of matching pursuit as target)
  VectorType Min(Dimension, 1e20);
  VectorType Max(Dimension, -1e20);

  for (int i = 0; i < m_NumCells; i++)
  {
    VectorType p = m_Centers.get_row(i);
    for (unsigned int dim=0; dim<Dimension; dim++)
    {
      Min[dim] = (p[dim] < Min[dim])?p[dim]:Min[dim];
      Max[dim] = (p[dim] > Max[dim])?p[dim]:Max[dim];
    }
  }

  Superclass::Superclass::m_BoundingBox.set_column(0, Min);
  Superclass::Superclass::m_BoundingBox.set_column(1, Max);
}

template class OrientedVolumeMesh<double,2>;
template class OrientedVolumeMesh<double,3>;
