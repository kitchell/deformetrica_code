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

#include <vector>
#include <memory>
#include "AnatomicalCoordinateSystem.h"
#include "LinearAlgebra.h"

using namespace def::algebra;

/**
 *  \brief      Abstract Geometry for the definition of geometrical supports.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The AbstractGeometry class encodes an object embedded in the current 2D or 3D space. The child classes contain
 *              the metric appearing in the fidelity-to-data term.\n \n
 *              See AbstractGeometry::AbstractGeometryType for the list of available object types.
 */
template<class ScalarType, unsigned int Dimension>
class AbstractGeometry {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///	Possible type of deformable object.
  typedef enum {
    null,                      /// Null value.
    ParametricImage,           /// Parametric image with gaussian kernel (see ParametricImage).
    SSDImage,                  /// Image with linear interpolation of voxel values and Sum Of Squared Differences (SSD) Metric
    LCCImage,                  /// Image with linear interpolation of voxel values and Local Correlation Coefficient (LCC) Metric
    EQLAImage,                 /// Image with linear interpolation of voxel values and "Ecart Quadratique au modele Local Affine" (EQLA) Metric (variant of LCC metric)
    MutualInformationImage,    /// Image with linear interpolation of voxel values and Mutual Information Metric (CAUTION : no gradient)
    Landmark,                  /// Landmark (see Landmark).
    OrientedPolyLine,          /// Current representation of a curve (see OrientedPolyLine).
    OrientedSurfaceMesh,       /// Current representation of a surface (see OrientedSurfaceMesh).
    NonOrientedPolyLine,       /// Varifold representation of a curve (see NonOrientedPolyLine).
    NonOrientedSurfaceMesh,    /// Varifold representation of a surface (see NonOrientedSurfaceMesh).
    PointCloud,                /// Point cloud (see PointCloud).
    OrientedVolumeMesh
  } AbstractGeometryType;

  /// Anatomical Coordinate System type.
  typedef AnatomicalCoordinateSystem<ScalarType, Dimension> AnatomicalCoordinateSystemType;



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  AbstractGeometry();
  /// Copy constructor.
  AbstractGeometry(const AbstractGeometry &other);
  /// Copy constructor with an update of the landmark point coordinates or a resampling of the image
  /// according to the new voxel positions in \e LandmarkOrImagePoints
  AbstractGeometry(const AbstractGeometry &example,
                   const MatrixType &LandmarkOrImagePoints);

  /// Returns a deformed version of the object
  std::shared_ptr<AbstractGeometry> DeformedObject(const MatrixType &LandmarkOrImagePoints) const {
    return doDeformedObject(LandmarkOrImagePoints);
  }
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<AbstractGeometry> doDeformedObject(const MatrixType &LandmarkOrImagePoints) const = 0;
 public:

  /// Makes a copy of the object.
  std::shared_ptr<AbstractGeometry> Clone() const { return doClone(); };
  /// Auxiliary private virtual function.
 private:
  virtual std::shared_ptr<AbstractGeometry> doClone() const = 0;
 public:

  /// Destructor.
  virtual ~AbstractGeometry();



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

// ScalarType GetTimePoint() { return m_Timepoint; }
// void SetTimePoint(ScalarType t) { m_Timepoint = t; }

// unsigned int GetTimeIndex() { return m_TimeIndex; }
// void SetTimeIndex(unsigned int t) { m_TimeIndex = t; }

  /// Returns true if the parameters of the deformation have changed, false otherwise.
  bool IsModified() const { return m_Modified; }
  /// Sets m_Modified to true (i.e. trajectory not computed or parameters changed).
  void SetModified() { m_Modified = true; }
  /// Sets m_Modified to false (i.e. trajectory computed and parameters not changed).
  void UnSetModified() { m_Modified = false; }

  /// Return the bounding box of the data
  MatrixType GetBoundingBox() const { return m_BoundingBox; }

  /// Returns the type of the deformable object.
  inline AbstractGeometryType GetType() const { return m_Type; }
  /// Returns true if the deformable object is of image kind, false otherwise.
  inline bool IsOfImageKind() const {
    return ((m_Type == SSDImage) || (m_Type == LCCImage) ||
        (m_Type == EQLAImage) || (m_Type == ParametricImage) ||
        (m_Type == MutualInformationImage));
  }
  /// Returns true if the deformable object is of landmark kind, false otherwise.
  inline bool IsOfLandmarkKind() const {
    return ((m_Type == Landmark) || (m_Type == PointCloud) ||
        (m_Type == OrientedPolyLine) || (m_Type == OrientedSurfaceMesh) ||
        (m_Type == NonOrientedPolyLine) || (m_Type == NonOrientedSurfaceMesh) ||
        (m_Type == OrientedVolumeMesh));
  }
  /// Returns true if the deformable object is a mesh, false otherwise.
  inline bool IsMesh() const {
    return ((m_Type == OrientedPolyLine) || (m_Type == NonOrientedPolyLine) ||
        (m_Type == OrientedSurfaceMesh) || (m_Type == NonOrientedSurfaceMesh) ||
        (m_Type == OrientedVolumeMesh));
  }
  /// Returns true if the deformable object is a parametric image, false otherwise.
  bool IsParametricImage() const { return (m_Type == ParametricImage); }
  /// Returns true if the deformable object is a linear interpolation image, false otherwise.
  bool IsLinearInterpImage() const {
    return ((m_Type == SSDImage) || (m_Type == LCCImage) ||
        (m_Type == EQLAImage));
  }
  /// Sets the type of the deformable object to ParametricImage.
  void SetParametricImageType() { m_Type = ParametricImage; }
  /// Sets the type of the deformable object to SSDImage.
  void SetSSDImageType() { m_Type = SSDImage; }
  /// Sets the type of the deformable object to LCCImage.
  void SetLCCImageType() { m_Type = LCCImage; }
  /// Sets the type of the deformable object to EQLAImage.
  void SetEQLAImageType() { m_Type = EQLAImage; }
  /// Sets the type of the deformable object to MutualInformationImage.
  void SetMutualInformationImageType() { m_Type = EQLAImage; }
  /// Sets the type of the deformable object to Landmark.
  void SetLandmarkType() { m_Type = Landmark; }
  /// Sets the type of the deformable object to PointCloud.
  void SetPointCloudType() { m_Type = PointCloud; }
  /// Sets the type of the deformable object to OrientedPolyLine.
  void SetOrientedPolyLineType() { m_Type = OrientedPolyLine; }
  /// Sets the type of the deformable object to OrientedSurfaceMesh.
  void SetOrientedSurfaceMeshType() { m_Type = OrientedSurfaceMesh; }
  /// Sets the type of the deformable object to NonOrientedPolyLine.
  void SetNonOrientedPolyLineType() { m_Type = NonOrientedPolyLine; }
  /// Sets the type of the deformable object to NonOrientedSurfaceMesh.
  void SetNonOrientedSurfaceMeshType() { m_Type = NonOrientedSurfaceMesh; }
  //Sets the type of the deformable object to OrientedVolumeMesh.
  void SetOrientedVolumeMeshType() { m_Type = OrientedVolumeMesh; }

  /// Returns the number of points of the object.
  virtual unsigned long GetNumberOfPoints() const = 0;

  /// Returns the dimension of the discretization of the object, used in Bayesian atlas update (for noise dimension).
  virtual unsigned long GetDimensionOfDiscretizedObject() const = 0;

  /// Returns the label of the anatomical orientation.
  std::string GetAnatomicalCoordinateSystemLabel() const { return m_AnatomicalOrientation.GetAnatomicalCoordinateSystemLabel(); }
  /// Sets the label of the anatomical orientation to \e label.
  void SetAnatomicalCoordinateSystem(std::string label) {
    m_AnatomicalOrientation.SetAnatomicalCoordinateSystemLabel(label);
  }
  /// Returns the change of basis matrix of the anatomical orientation.
  MatrixType GetAnatomicalCoordinateSystemMatrix() const { return m_AnatomicalOrientation.GetChangeOfBasisMatrix(); }
  /// Sets the change of basis matrix of the anatomical orientation to \e m.
  void SetAnatomicalCoordinateSystem(MatrixType m) { m_AnatomicalOrientation.SetChangeOfBasisMatrix(m); }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Updates several information concerning the deformable object (see details on child class).
  virtual void Update() = 0;

  /**
   *  \brief      Returns the norm between itself and the target.
   *
   *  \details    Denoting by \e S the source (i.e. the class itself) and by \e T the target object, this
   *              method computes the following norm :
   *              \f[
   *              \left\Vert S - T \right\Vert_{W}^2 ,
   *              \f]
   *              where the metric \f$ \left\Vert ~.~ \right\Vert_{W}^2 \f$ depends on the child class.
   *
   *  \param[in]  target	Target of same type as the deformable object.
   *  \return     The norm of the difference.
   */
  virtual ScalarType ComputeMatch(const std::shared_ptr<AbstractGeometry> target) = 0;

  /**
   *  \brief      Returns the gradient of the norm between itself and the target.
   *
   *  \details    Denoting by \e S the source (i.e. the class itself) and by \e T the target object, this
   *              method computes the following gradient :
   *              \f[
   *              \nabla\left\Vert S - T \right\Vert_{W}^2 ,
   *              \f]
   *              where the metric \f$ \left\Vert ~.~ \right\Vert_{W}^2 \f$ depends on the child class.
   *
   *  \param[in]  target	Target of same type as the deformable object.
   *  \return	    The gradient of the norm of the difference.
   */
  virtual MatrixType ComputeMatchGradient(const std::shared_ptr<AbstractGeometry> target) = 0;

  /// Saves in \e filename the deformable object (for Landmark type, it is saved in *.vtk format (VTK PolyData)).
  virtual void WriteObject(std::string filename) const = 0;
  virtual void WriteObject(std::string filename, const MatrixListType &velocity) const {
    WriteObject(filename);
  } // Further outputs could be produced. TODO.


 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Computes the union of two bounding box.
  MatrixType UnionBoundingBox(MatrixType BB1, MatrixType BB2);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Type of the deformable object.
  AbstractGeometryType m_Type;

  /// Boolean which avoids re-updating computation of centers, tangents and self-norm
  bool m_Modified;

  /// Tight Bounding box, which contain the data (i.e. centers of mesh cells and centers of voxels)
  MatrixType m_BoundingBox;

  /// Anatomical coordinate system associated to the deformable object.
  AnatomicalCoordinateSystemType m_AnatomicalOrientation;

}; /* class AbstractGeometry */

