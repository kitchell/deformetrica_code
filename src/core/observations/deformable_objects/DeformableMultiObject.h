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

/// Core files.
#include "AbstractGeometry.h"
#include "Landmark.h"
#include "LinearInterpImage.h"
#include "ParametricImage.h"

/// Support files.
#include "LinearAlgebra.h"
#include "KernelFactory.h"
#include "GridFunctions.h"

/// Librairies files.
#include <vector>
#include "itkImage.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"

using namespace def::algebra;

/**
 *  \brief      Collections of deformable objects (i.e. multi-object)
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The DeformableMultiObject class is used to deal with collections of deformable objects
 *              embedded in the current 2D or 3D space. It extends for such collections the methods
 *              of the deformable object class that are designed for a single object at a time.
 */
template <class ScalarType, unsigned int Dimension>
class DeformableMultiObject {
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Deformable object type.
	typedef AbstractGeometry<ScalarType, Dimension> AbstractGeometryType;

	/// List of deformable objects type.
	typedef std::vector<std::shared_ptr<AbstractGeometryType>> AbstractGeometryList;

	/// ITK image type.
	typedef itk::Image<ScalarType, Dimension> ImageType;
	/// ITK image pointer type.
	typedef typename ImageType::Pointer ImageTypePointer;

    /// Landmark type.
    typedef Landmark<ScalarType, Dimension> LandmarkType;
    /// Linear interpolation image type.
    typedef LinearInterpImage<ScalarType, Dimension> LIImageType;
    /// Parametric image type (gaussian kernel).
    typedef ParametricImage<ScalarType, Dimension> ParametricImageType;

    /// Grid functions type.
    typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Constructor.
	DeformableMultiObject();

	/// Copy constructor.
	DeformableMultiObject(const DeformableMultiObject& other);

	/// Contructor which copies the \e Example and sets the coordinates of landmark points to \e LandmarkPoints
	/// and re-sample to image according to voxels positions in \e ImagePoints. It is used to get
	/// the deformed version of a DeformableMultiObject.
	DeformableMultiObject(const DeformableMultiObject& Example,
						  const MatrixType& LandmarkPoints,
						  const MatrixType& DownSampledImageMap);
	
	/// Returns a deformed version of the object.
	std::shared_ptr<DeformableMultiObject> DeformedMultiObject(MatrixType& LandmarkPoints,
															   MatrixType& ImagePoints) const {
		return doDeformedMultiObject(LandmarkPoints, ImagePoints); }

	/// Makes a copy of the object.
	std::shared_ptr<DeformableMultiObject> Clone() const { return std::make_shared<DeformableMultiObject>(*this); }

	/// Destructor.
	~DeformableMultiObject();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Return the list of the deformable objects.
	inline AbstractGeometryList GetObjectList() const { return m_ObjectList; }

	/// Set the list of deformable objects as a copy of \e list.
	void SetObjectList(AbstractGeometryList& list);

	/// Returns a vector of bools, whose ith element is true iff m_ObjectList[k] is of image kind.
	std::vector<bool> IsOfImageKind() const;
	/// Returns a vector of bools, whose ith element is true iff m_ObjectList[k] is of landmark kind.
	std::vector<bool> IsOfLandmarkKind() const;

	/// Returns the coordinates of the vertices for Landmark type (or children object).
	MatrixType GetLandmarkPoints() const;
	/// Returns the voxel positions of the downsample image in LinearInterpImage type objects (or children object).
	MatrixType GetDownSampledImageMap() const;
	/// Returns an example of downsampled image, whose voxel coordinates are given in this->GetDownSampledImageMap()
	ImageTypePointer GetDownSampledImage() const;
	/// Returns image at full resolution
	ImageTypePointer GetImage() const;
	
	// /// Set new coordinates of all landmark type objects (used mostly to compute a deformed object)
	// void SetLandmarkPoints(const MatrixType& Y);
	// /// Set new positions of the voxel in the image domain (used mostly to compute a deformed image)
	// void SetImagePoints(const MatrixType& Y);

	/// Return the coordinates of the vertices of the landmark-type objects and the intensities of the image object, if any (used mostly for updating template parameters)
	MatrixListType GetImageIntensityAndLandmarkPointCoordinates() const;
	/// Set new vertex coordinates in landmark type objects and new image intensities in the image type object (used mostly for updating template parameters)
	void UpdateImageIntensityAndLandmarkPointCoordinates(const MatrixListType & Y);

	/// Return the index of the image in the list of objects (returns -1 if there is no image).
	int GetImageIndex() const { return m_ImageIndex; }

	/// Return the number of objects of image kind.
	int GetNumberOfImageKindObjects() const { return m_NumberOfImageKindObjects; }
	/// Return the number of objects of landmark kind.
	int GetNumberOfLandmarkKindObjects() const { return m_NumberOfLandmarkKindObjects; } 
	/// Return the total number of objects.
	int GetNumberOfObjects() const { return m_NumberOfObjects; }

	/// Returns the list of each number of points associated to the deformable objects.
	std::vector<int> GetNumberOfPoints() const { return m_NumberOfPoints; }
	/// Returns the number of points of all deformable objects.
	int GetTotalNumberOfPoints() const { return (m_NumberOfLandmarkPoints + m_NumberOfImagePoints); }

	/// Returns the number of landmark points
	int GetNumberOfLandmarkPoints() const { return m_NumberOfLandmarkPoints; }
	
	/// Returns the number of image points
	int GetNumberOfImagePoints() const { return m_NumberOfImagePoints; }
	
	/// Returns a tight bounding box enclosing all objects
	MatrixType GetBoundingBox() const { return m_BoundingBox; };
	/// Computes a tight bounding box enclosing all objects.
	void UpdateBoundingBox();
	
	/// Get dimension of discretized objects for each object
	std::vector<unsigned long> GetDimensionOfDiscretizedObjects() const;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates information about the deformable objects and updates the bounding box.
	void Update();

	/// Computes AbstractGeometry::ComputeMatch(AbstractGeometry* target) for each deformable object.
	std::vector<ScalarType> ComputeMatch(const std::shared_ptr<DeformableMultiObject> target);

	/// Computes AbstractGeometry::ComputeMatchGradient(AbstractGeometry* target) for each deformable object.
	MatrixListType ComputeMatchGradient(const std::shared_ptr<DeformableMultiObject> target);

	/// Transforms concatenated data located at landmark points and image points into a list
	void ListToMatrices(const MatrixListType& L, MatrixType& MLandmark, MatrixType& MImage) const;
	/// Transforms list of objects into concatenated data of landmark types and image types
	void MatricesToList(MatrixListType& L, const MatrixType& MLandmark, const MatrixType& MImage) const;
	
	/// Splats the difference between image in this and image in \e targ
	MatrixType SplatDifferenceImage(const std::shared_ptr<DeformableMultiObject> targ,
                                    const MatrixType& downSampledImageMap) const;

	/// Calls the method AbstractGeometry::WriteObject() for each deformable object.
	void WriteMultiObject(std::vector<std::string>& str) const;
	void WriteMultiObject(std::vector<std::string>& str, const MatrixListType& velocity) const;


protected:

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Protected method(s).
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Returns a deformed version of the object.
    virtual std::shared_ptr<DeformableMultiObject> doDeformedMultiObject(MatrixType& LandmarkPoints,
                                                                         MatrixType& ImagePoints) const {
        return std::make_shared<DeformableMultiObject>(*this, LandmarkPoints, ImagePoints); }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Collection of m_NumberOfObjects deformable objects.
	AbstractGeometryList m_ObjectList;

	/// Total number of objects.
	int m_NumberOfObjects;

	/// Number of objects of image kind.
	int m_NumberOfImageKindObjects;

	/// Number of objects of landmark kind.
	int m_NumberOfLandmarkKindObjects;

	/// Vector containing at each cell the number of points of the associated deformable object.
	std::vector<int> m_NumberOfPoints;
	
	/// Total number of vertices in the multi-object
	int m_NumberOfLandmarkPoints;
	
	/// Number of pixels/voxels in the image
	int m_NumberOfImagePoints;
	
	/// Index of the object of image type (at most one image allowed within a multi-object).
	int m_ImageIndex;
	
	/// Bounding Box tightly enclosing all objects
	MatrixType m_BoundingBox;


}; /* class DeformableMultiObject */
