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

#include "LinearAlgebra.h"
#include "KernelFactory.h"
#include <vector>

#include "AbstractGeometry.h"
#include "DeformableMultiObject.h"

#include "itkImage.h"

#include <memory>

using namespace def::algebra;

/**
 *  \brief      Abstract deformations.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The AbstractDeformations class encodes diffeomorphic deformations of the ambient 2D or 3D space,
 *              which then deforms objects embedded into it.
 */
template<class ScalarType, unsigned int Dimension>
class AbstractDeformations {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///	Possible type of deformations.
  typedef enum {
    null,            /*!< Default value. */
    Diffeos            /*!< Standard deformation (see Diffeos). */
  } DeformationsType;

  /// Deformable object type.
  typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;
  /// Vector of deformable objects type.
  typedef std::vector<AbstractGeometryType *> AbstractGeometryList;
  /// Deformable multi-object type.
  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;
  /// Template type.
  // typedef Template<ScalarType, Dimension> TemplateType;

  /// ITK image type.
  typedef itk::Image<ScalarType, Dimension> ImageType;
  /// ITK image pointer type.
  typedef typename ImageType::Pointer ImageTypePointer;



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Contructor.
  AbstractDeformations();

  ///	Copy constructor.
  AbstractDeformations(const AbstractDeformations &other);

  ///	Makes a copy of the object.
  std::shared_ptr<AbstractDeformations> Clone() { return doClone(); }

  /// Destructor.
  virtual ~AbstractDeformations();

  // // copy essential information, i.e. everything but the deformable multi-object
  // virtual void CopyInformation(const AbstractDeformations& other);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns true if it is of Diffeos type, false otherwise.
  inline bool isADiffeo() const { return (m_Type == Diffeos); }
  /// Sets type of the deformation to Diffeos.
  inline void SetDiffeosType() { m_Type = Diffeos; }

  /// Returns true if one of the parameters of the deformation has changed, false otherwise.
  inline bool IsModified() const { return m_Modified; }
  /// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
  inline void SetModified() { m_Modified = true; }
  /// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
  inline void UnsetModified() { m_Modified = false; }

  /// Gets the collection of deformable objects.
  std::shared_ptr<DeformableMultiObjectType> GetDeformableMultiObject() const { return m_DeformableMultiObject; }
  /// Sets the collection of deformable objects to \e DMO.
  void SetDeformableMultiObject(std::shared_ptr<DeformableMultiObjectType> DMO);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Computes the deformation of the deformable objects under the action of the child class of AbstractDeformations.
  virtual void Update() = 0;

  /// Returns the list of deformed objects at final time.
  virtual std::shared_ptr<DeformableMultiObjectType> GetDeformedObject() const = 0;

  /// Saves the flow of each deformable object.
  virtual void WriteFlow(const std::vector<std::string> &name, const std::vector<std::string> &extension) = 0;

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///	Makes a copy of the object.
  virtual std::shared_ptr<AbstractDeformations> doClone() = 0;

  /// Computes trajectories of voxel positions using inverse deformation (used to map the image domain with \f$\phi^{-1}\f$).
  virtual void FlowImagePointsTrajectory() = 0;
  /// Computes trajectories of any points in space (used for all objects of Landmark type or children types).
  virtual void FlowLandmarkPointsTrajectory() = 0;



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Type of the deformation.
  DeformationsType m_Type;

  /// Pointer on the collection of deformable objects on which our deformation will act.
  std::shared_ptr<DeformableMultiObjectType> m_DeformableMultiObject;

  /// Set of points, which are transformed by the deformation \f$\phi(x)\f$
  MatrixType m_LandmarkPoints;
  /// Set of points, which are transformed by the inverse deformation \f$\phi^{-1}(x)\f$ (typically the voxel's positions in \e m_DownSampledWorkingImage)
  MatrixType m_ImagePoints;
  /// Downsampled reference image used to deform the image domain
  ImageTypePointer m_DownSampledImage;
  /// Full resolution image used essentially in adjoint equations
  ImageTypePointer m_Image;

  /// true if there is at least one object of landmark type (or children types) in the multi-object
  bool m_IsLandmarkPoints;
  /// true if there is an image in the multi-object
  bool m_IsImagePoints;

  /// Boolean which requires to update deformation parameters
  bool m_Modified;
  ///	Boolean which avoids computing the trajectory (via Update()) if no parameter has changed.
  bool m_DeformableObjectModified;

}; /* class AbstractDeformations */

