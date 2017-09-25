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

#include "AbstractGeometry.h"

#include "itkSimpleFastMutexLock.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkVersion.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief 		Landmarks (i.e. labelled point sets)
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The Landmark class inherited from AbstractGeometry represents a set of labelled points.
 *              This class assumes that the source and the target have the same number of points with a
 *              point-to-point correspondence.
 */
template <class ScalarType, unsigned int Dimension>
class Landmark : public AbstractGeometry<ScalarType, Dimension>
{
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Deformable object type.
  typedef AbstractGeometry<ScalarType, Dimension> Superclass;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  Landmark();

  /// Copy constructor.
  Landmark(const Landmark& other);
  /// Constructor which copies the objects and update vertex coordinates.
  Landmark(const Landmark& example, const MatrixType& LandmarkPoints);

  /// Returns a deformed version of the point set, where the deformation is given by the position of vertices in \e LandmarkPoints.
  std::shared_ptr<Landmark> DeformedObject(const MatrixType& LandmarkPoints) const {
    return std::static_pointer_cast<Landmark>(doDeformedObject(LandmarkPoints)); }
 private:
  virtual std::shared_ptr<Superclass> doDeformedObject(const MatrixType& LandmarkPoints) const {
    return std::static_pointer_cast<Superclass>(std::make_shared<Landmark>(*this, LandmarkPoints)); }
 public:

  /// Clones the object.
  std::shared_ptr<Landmark> Clone() const { return std::static_pointer_cast<Landmark>(doClone()); }
 private:
  virtual std::shared_ptr<Superclass> doClone() const {
    return std::static_pointer_cast<Superclass>(std::make_shared<Landmark>(*this)); }
 public:

  /// Destructor.
  virtual ~Landmark();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the pointer on a VTK object to a reference \e polyData.
  virtual void SetPolyData(vtkPointSet* polyData);
  /// Transforms polyData from native orientation system into ITK (LPS) default orientation,
  /// or inversely, if the flag "inverse" is set to true.
  vtkSmartPointer<vtkPointSet> ReorientPolyData(vtkPointSet* polyData, bool inverse) const;

  /// Update the PolyData with new coordinates of vertices. Need a call to Update() afterwards.
  void UpdatePolyDataWithPointCoordinates(const MatrixType& LandmarkPoints);

  /// Returns the vertex coordinates as a matrix.
  MatrixType GetPointCoordinates() const { return m_PointCoordinates; }

  /// Returns the number of points of the deformable object.
  virtual unsigned long GetNumberOfPoints() const { return m_NumberOfPoints; }

  /// Returns the dimension of the discretized image, here the number of points times the dimension.
  virtual unsigned long GetDimensionOfDiscretizedObject() const { return (Dimension * m_NumberOfPoints); }
  /// Given a point, this method tries to identify the cell it belongs to. Throws an error if it fails. The method also gets the weights for interpolation.
  virtual int GetCellForPosition(const VectorType& pos, double* weights);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Update();

  virtual ScalarType ComputeMatch(const std::shared_ptr<Superclass> target);

  virtual MatrixType ComputeMatchGradient(const std::shared_ptr<Superclass> target);

  virtual void WriteObject(std::string filename) const;

  ///True if the data used if of type Unstructured Grid (OrientedVolumeMesh for now)
  inline bool IsOfUnstructuredKind() const {return (this->m_Type == Superclass::AbstractGeometryType::OrientedVolumeMesh) && Dimension == 3;}


 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates the bounding box of the data.
  void UpdateBoundingBox();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///	Pointer on a geometric structure namely a VTK object.
  vtkSmartPointer<vtkPointSet> m_PointSet;

  ///	Matrix coordinates of the points (Size : NumberOfPoints x Dimension).
  MatrixType m_PointCoordinates;

  ///	Number of points of the deformable object.
  int m_NumberOfPoints;

  ///	Object used to perform mutex (mutual exclusion) with m_ReferencePolyData (important for multithreaded programming).
  itk::SimpleFastMutexLock m_VTKMutex;


}; /* class Landmark */

