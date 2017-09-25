/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "OrientedSurfaceMesh.h"

#include "KernelFactory.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkTriangleFilter.h"
#include "myvtkPolyDataNormals.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include <cstring>
#include <iostream>
#include <sstream>

#include <cassert>



////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
OrientedSurfaceMesh<ScalarType, Dimension>
::OrientedSurfaceMesh() : Superclass()
{
  if (Dimension != 3)
    throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

  this->SetOrientedSurfaceMeshType();
  m_KernelWidth = 0;
  m_KernelType = null;
  m_Reorient = false;
}



template <class ScalarType, unsigned int Dimension>
OrientedSurfaceMesh<ScalarType, Dimension>
::~OrientedSurfaceMesh()
{}



template <class ScalarType, unsigned int Dimension>
OrientedSurfaceMesh<ScalarType, Dimension>
::OrientedSurfaceMesh(const OrientedSurfaceMesh& other) : Superclass(other)
{
  if (Dimension != 3)
    throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

  this->SetOrientedSurfaceMeshType();

  m_Centers = other.m_Centers;
  m_Normals = other.m_Normals;
  m_NumCells = other.m_NumCells;

  m_KernelWidth = other.m_KernelWidth;
  m_KernelType = other.m_KernelType;
  m_NormSquared = other.m_NormSquared;

  m_Reorient = other.m_Reorient;
}



template <class ScalarType, unsigned int Dimension>
OrientedSurfaceMesh<ScalarType, Dimension>
::OrientedSurfaceMesh(const OrientedSurfaceMesh& ex, const MatrixType& LP) : Superclass(ex, LP)
{
  if (Dimension != 3)
    throw std::runtime_error("Cannot create Surface Mesh Object in Dimension other than 3");

  this->SetOrientedSurfaceMeshType();

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
OrientedSurfaceMesh<ScalarType, Dimension>
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
ScalarType OrientedSurfaceMesh<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<AbstractGeometryType> target)
{

  if (this->GetType() != target->GetType())
    throw std::runtime_error("AbstractGeometries types mismatched");

  const std::shared_ptr<const OrientedSurfaceMesh> targetOrientedSurfaceMesh
      = std::static_pointer_cast<const OrientedSurfaceMesh>(target);

  if (m_KernelWidth != targetOrientedSurfaceMesh->GetKernelWidth())
    throw std::runtime_error("Kernel width of surface currents mismatched");

  this->Update();
  // target->Update();

  MatrixType targCenters = targetOrientedSurfaceMesh->GetCenters();
  MatrixType targNormals = targetOrientedSurfaceMesh->GetNormals();

  ScalarType match = targetOrientedSurfaceMesh->GetNormSquared() + this->GetNormSquared();

  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

  MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

  std::shared_ptr<KernelType> kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
  kernelObject->SetKernelWidth(m_KernelWidth);
  kernelObject->SetSources(m_Centers);
  kernelObject->SetWeights(m_Normals);

  MatrixType SdotT = kernelObject->Convolve(targCenters);
  for (int i = 0; i < targetOrientedSurfaceMesh->GetNumberOfCells(); i++)
    match -= 2.0f * dot_product(SdotT.get_row(i), targNormals.get_row(i));

  return match;
}


template<class ScalarType, unsigned int Dimension>
MatrixType OrientedSurfaceMesh<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target)
{
  if (this->GetType() != target->GetType())
    throw std::runtime_error("Deformable objects types mismatched");

  const std::shared_ptr<const OrientedSurfaceMesh> targetOrientedSurfaceMesh
      = std::static_pointer_cast<const OrientedSurfaceMesh>(target);

  if (m_KernelWidth != targetOrientedSurfaceMesh->GetKernelWidth())
    throw std::runtime_error("Kernel width of surface currents mismatched");

  this->Update();
  // target->Update();

  MatrixType targCenters = targetOrientedSurfaceMesh->GetCenters();
  MatrixType targNormals = targetOrientedSurfaceMesh->GetNormals();

  MatrixType Pts = this->GetPointCoordinates();

  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

  MatrixType DataDomain = this->UnionBoundingBox(this->GetBoundingBox(), target->GetBoundingBox());

  std::shared_ptr<KernelType> kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), DataDomain);
  kernelObject->SetKernelWidth(m_KernelWidth);

  kernelObject->SetSources(m_Centers);
  kernelObject->SetWeights(m_Normals);
  MatrixType KtauS = kernelObject->Convolve(m_Centers);
  MatrixType gradKtauS = kernelObject->ConvolveGradient(m_Centers, m_Normals);

  kernelObject->SetSources(targCenters);
  kernelObject->SetWeights(targNormals);
  MatrixType KtauT = kernelObject->Convolve(m_Centers);
  MatrixType gradKtauT = kernelObject->ConvolveGradient(m_Centers, m_Normals);

  MatrixType gradKtau = (gradKtauS - gradKtauT) * 2.0f / 3.0f;

  MatrixType gradmatch(this->GetNumberOfPoints(), Dimension, 0.0);
  for (int f = 0; f < m_NumCells; f++)
  {
    vtkSmartPointer < vtkIdList > ptIds = vtkSmartPointer<vtkIdList>::New();

    m_VTKMutex.Lock();
    Superclass::m_PointSet->GetCellPoints(f, ptIds);
    m_VTKMutex.Unlock();

    int ind0 = ptIds->GetId(0);
    int ind1 = ptIds->GetId(1);
    int ind2 = ptIds->GetId(2);

    VectorType e0 = Pts.get_row(ind2) - Pts.get_row(ind1);
    VectorType e1 = Pts.get_row(ind0) - Pts.get_row(ind2);
    VectorType e2 = Pts.get_row(ind1) - Pts.get_row(ind0);

    VectorType Ktau = KtauS.get_row(f) - KtauT.get_row(f);
    gradmatch.set_row(ind0, gradmatch.get_row(ind0) + cross_3d(e0, Ktau));
    gradmatch.set_row(ind1, gradmatch.get_row(ind1) + cross_3d(e1, Ktau));
    gradmatch.set_row(ind2, gradmatch.get_row(ind2) + cross_3d(e2, Ktau));

    gradmatch.set_row(ind0, gradmatch.get_row(ind0) + gradKtau.get_row(f));
    gradmatch.set_row(ind1, gradmatch.get_row(ind1) + gradKtau.get_row(f));
    gradmatch.set_row(ind2, gradmatch.get_row(ind2) + gradKtau.get_row(f));

  }

  return gradmatch;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
OrientedSurfaceMesh<ScalarType, Dimension>
::CheckMeshAndNormals()
{
  //We downcast the pointer to a vtkPolyData first.
  vtkPolyData* polydata = dynamic_cast<vtkPolyData*>(Superclass::m_PointSet.GetPointer());

  vtkSmartPointer<vtkTriangleFilter> trif =
      vtkSmartPointer<vtkTriangleFilter>::New();
#if (VTK_MAJOR_VERSION == 5)
  trif->SetInput(Superclass::m_PolyData);
#else
  trif->SetInputData(polydata);
#endif
  trif->PassVertsOff();
  trif->PassLinesOff();
  trif->Update();

  vtkSmartPointer<myvtkPolyDataNormals> normalf = vtkSmartPointer<myvtkPolyDataNormals>::New();
#if (VTK_MAJOR_VERSION == 5)
  normalf->SetInput(trif->GetOutput());
#else
  normalf->SetInputData(trif->GetOutput());
#endif
  normalf->ComputePointNormalsOff();
  normalf->ComputeCellNormalsOn();
  normalf->SplittingOff();

  if (m_Reorient)
  {
    normalf->ConsistencyOn();
    normalf->AutoOrientNormalsOn(); // Should have closed surface
    normalf->FlipNormalsOn();
  }
  else
  {
    normalf->ConsistencyOff();
    normalf->AutoOrientNormalsOff();
  }
  normalf->Update();

  vtkSmartPointer<vtkPolyData> output = normalf->GetOutput();
  output -> BuildLinks();

  Superclass::m_PointSet = vtkSmartPointer<vtkPointSet>(output);

  m_NumCells = Superclass::m_PointSet->GetNumberOfCells();

}



template <class ScalarType, unsigned int Dimension>
void
OrientedSurfaceMesh<ScalarType, Dimension>
::UpdateCentersNormals()
{
  const MatrixType Pts = this->GetPointCoordinates();
  m_Centers.set_size(m_NumCells, Dimension);
  m_Normals.set_size(m_NumCells, Dimension);

  for (unsigned int i = 0; i < m_NumCells; i++)
  {
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    m_VTKMutex.Lock();
    Superclass::m_PointSet->GetCellPoints(i, ptIds);
    m_VTKMutex.Unlock();

    if (ptIds->GetNumberOfIds() != 3)
      throw std::runtime_error("Not a triangle cell!");

    VectorType p0 = Pts.get_row(ptIds->GetId(0));
    VectorType p1 = Pts.get_row(ptIds->GetId(1));
    VectorType p2 = Pts.get_row(ptIds->GetId(2));

    m_Centers.set_row(i, (p0 + p1 + p2) / 3.0f );
    m_Normals.set_row(i, cross_3d(p2 - p0, p1 - p0) / 2);

  }
}



template <class ScalarType, unsigned int Dimension>
void
OrientedSurfaceMesh<ScalarType, Dimension>
::UpdateSelfNorm()
{
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

  std::shared_ptr<KernelType> kernelObject = kFactory->CreateKernelObject(this->GetKernelType(), this->GetBoundingBox());

  kernelObject->SetKernelWidth(m_KernelWidth);
  kernelObject->SetSources(m_Centers);
  kernelObject->SetWeights(m_Normals);

  MatrixType selfKW = kernelObject->Convolve(m_Centers);

  m_NormSquared = 0;

  for (unsigned int i = 0; i < m_NumCells; i++)
    m_NormSquared += dot_product(selfKW.get_row(i), m_Normals.get_row(i));
}



template <class ScalarType, unsigned int Dimension>
void
OrientedSurfaceMesh<ScalarType, Dimension>
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

template class OrientedSurfaceMesh<double,2>;
template class OrientedSurfaceMesh<double,3>;

