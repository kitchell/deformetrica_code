/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

#include "Landmark.h"

#include "KernelType.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      Oriented surface meshes.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *  \details    The OrientedSurfaceMesh class inherited from Landmark represents a triangular mesh.
 *              This class uses the current representation of the surface, which is sensitive to the local orientation
 *              of the mesh (a flip in ordering of the vertex indices in the connectivity matrix changes the sign
 *              of the surface element (i.e. the triangle) in the space of currents).
 */
template <class ScalarType, unsigned int Dimension>
class OrientedSurfaceMesh : public Landmark<ScalarType, Dimension>
{
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Landmark type.
  typedef Landmark<ScalarType, Dimension> Superclass;

  /// Deformable object type.
  typedef typename Superclass::Superclass AbstractGeometryType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Contructor.
  OrientedSurfaceMesh();

  /// Copy constructor.
  OrientedSurfaceMesh(const OrientedSurfaceMesh& other);
  /// Constructor which copies the object and update vertex coordinates
  OrientedSurfaceMesh(const OrientedSurfaceMesh& example, const MatrixType& LandmarkPoints);

  /// Returns a Deformed version of the mesh, where the deformation is given by the position of vertices in \e LandmarkPoints
  std::shared_ptr<OrientedSurfaceMesh> DeformedObject(const MatrixType& LandmarkPoints) const {
      return std::static_pointer_cast<OrientedSurfaceMesh>(doDeformedObject(LandmarkPoints)); }

  /// Clones the object.
  std::shared_ptr<OrientedSurfaceMesh> Clone() const {
      return std::static_pointer_cast<OrientedSurfaceMesh>(doClone()); }

  /// Destructor.
  virtual ~OrientedSurfaceMesh();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Uses VTK filter to consistently re-orient normals for genus 0 surfaces.
  inline void SetReorient() { m_Reorient = true; this->SetModified(); }

  /// Unsets the use of VTK filter to automatically re-orient surface normals.
  inline void UnSetReorient() { m_Reorient = false; this->SetModified(); }

  /// Returns the centers of the cells.
  inline MatrixType GetCenters() const { return m_Centers; }

  /// Returns the normals of the cells.
  inline MatrixType GetNormals() const { return m_Normals; }

  /// Returns the number of cells.
  inline int GetNumberOfCells() const { return m_NumCells; }

  ///	Returns the type of the kernel.
  inline KernelEnumType GetKernelType() const { return m_KernelType; }
  /// Sets the type of the kernel to \e kernelType.
  inline void SetKernelType(KernelEnumType kernelType) { m_KernelType = kernelType; this->SetModified(); }

  ///	Returns the size of the kernel.
  inline ScalarType GetKernelWidth() const { return m_KernelWidth; }
  /// Sets the size of the kernel to \e h.
  inline void SetKernelWidth(ScalarType h) {	m_KernelWidth = h; this->SetModified(); }

  /// Returns the RKHS-norm of itself.
  inline ScalarType GetNormSquared() const { return m_NormSquared; }

  // STANLEY
  /*
  virtual int GetDimensionOfDiscretizedObject() const
  {
      int d = 0;
      for (unsigned int dim = 0; dim < Dimension; dim ++)
      {
          d += floor( (Superclass::Superclass::m_BoundingBox(dim,1) - Superclass::Superclass::m_BoundingBox(dim,0)) / m_KernelWidth );
      }
      d *= Dimension;
      return d;
  }
  */
  // PIETRO
  virtual unsigned long GetDimensionOfDiscretizedObject() const
  {
      int d = 1;
      for (unsigned int dim = 0; dim < Dimension; dim ++)
      {
          d *= floor( (Superclass::Superclass::m_BoundingBox(dim,1) - Superclass::Superclass::m_BoundingBox(dim,0)) / m_KernelWidth + 1.0);
      }
      d *= Dimension;
      return d;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Update();

  /// See AbstractGeometry::ComputeMatch(AbstractGeometry* target) for details.
  virtual ScalarType ComputeMatch(const std::shared_ptr<AbstractGeometryType> target);

  /// See AbstractGeometry::ComputeMatchGradient(AbstractGeometry* target) for details.
  virtual MatrixType ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target);



 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns a Deformed version of the mesh, where the deformation is given by the position of vertices in \e LandmarkPoints
  virtual std::shared_ptr<AbstractGeometryType> doDeformedObject(const MatrixType& LandmarkPoints) const {
      return std::static_pointer_cast<AbstractGeometryType>(std::make_shared<OrientedSurfaceMesh>(*this, LandmarkPoints)); }

  /// Clones the object.
  virtual std::shared_ptr<AbstractGeometryType> doClone() const {
      return std::static_pointer_cast<AbstractGeometryType>(std::make_shared<OrientedSurfaceMesh>(*this)); }


  /// Updates the bounding box.
  void UpdateBoundingBox();

  /// Possibly reorient normals according to vtk filters.
  /// \warning   Make sure the ordering of point cells is consistent with the direction of the normal.
  void CheckMeshAndNormals();

  /// Updates the centers and the normals from the points.
  void UpdateCentersNormals();

  /*
   *	\brief		Computes the centers and the normals from the points.
   *
   *	\details	Given a triangular mesh, this method computes from the vertices of the mesh
   *				the centers and the normals of each face.
   *
   *	\param[in]	Pts			The vertices of the mesh.
   *	\param[out]	Centers		The centers of the cells (Size : NumCells x Dimension).
   *	\param[out]	Normals		The normals of the cells (Size : NumCells x Dimension).
   */
  //void ComputeCentersNormals(const MatrixType& Pts, MatrixType& Centers, MatrixType& Normals);

  /// Computes the RKHS-norm of itself.
  void UpdateSelfNorm();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///	Matrix coordinates of the centers of the cells  (Size : NumCells x Dimension).
  MatrixType m_Centers;

  ///	Matrix coordinates of the normals of the cells  (Size : NumCells x Dimension).
  MatrixType m_Normals;

  /// Number of cells (i.e. triangles) of the mesh.
  int m_NumCells;

  /// true to use VTK filter to re-orient normals of genus-0 surfaces
  bool m_Reorient;

  ///	Type of the kernel.
  KernelEnumType m_KernelType;
  ///	Size of the kernel.
  ScalarType m_KernelWidth;

  /// Squared RKHS-norm of the oriented surface.
  ScalarType m_NormSquared;

  /// See Landmark::m_VTKMutex for details.
  itk::SimpleFastMutexLock m_VTKMutex;


}; /* class OrientedSurfaceMesh */
