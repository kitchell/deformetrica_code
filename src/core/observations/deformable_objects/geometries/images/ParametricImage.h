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

/// Class file.
#include "AbstractGeometry.h"

/// Core files.
#include "LinearInterpImage.h"

/// Support files.
#include "KernelFactory.h"
#include "GridFunctions.h"

/// Librairies files.
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>


/**
 *  \brief      An parametric image using kernel interpolation of control points values.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The ParametricImage class inherited from DeformableObject models continuous images
 *              (i.e. with infinite resolution) by interpolating control point intensities through a kernel.
 */
template<class ScalarType, unsigned int Dimension>
class ParametricImage : public AbstractGeometry<ScalarType, Dimension>
{
public:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// typedef 
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Deformable Object type.
	typedef AbstractGeometry<ScalarType, Dimension> Superclass;

    /// ITK image type.
    typedef itk::Image<ScalarType, Dimension> ImageType;
    /// ITK image pointer type.
    typedef typename ImageType::Pointer ImagePointerType;
    /// ITK duplicator type.
    typedef typename itk::ImageDuplicator<ImageType> DuplicatorType;
    /// ITK iterator type.
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    /// List of ITK images.
    typedef std::vector<ImagePointerType> ImagePointerListType;

    /// Kernel factory type.
    typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
    /// Kernel type.
    typedef typename KernelFactoryType::KernelBaseType KernelType;

	/// Linear interpolation image type.
	typedef LinearInterpImage<ScalarType, Dimension> LIImageType;
	/// Grid functions type.
	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor(s) / Destructor :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Constructor.
	ParametricImage();

	/// Copy constructor.
	ParametricImage(const ParametricImage& other);
    /// Copy constructor which modifies the geometry (encoded by new voxel positions \e imagePoints).
	ParametricImage(const ParametricImage& example, const MatrixType& imagePoints);

    /// Returns a deformed version of the image, where the deformation is given by the position of voxels positions in \e imagePoints.
	std::shared_ptr<ParametricImage> DeformedObject(const MatrixType& imagePoints) const {
		return std::static_pointer_cast<ParametricImage>(doDeformedObject(imagePoints)); }

	/// Makes a copy of the object.
	std::shared_ptr<ParametricImage> Clone() const { return std::static_pointer_cast<ParametricImage>(doClone()); }

	/// Destructor.
	virtual ~ParametricImage();


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Encapsulation method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Gets the photometric control points.
	MatrixType GetPhotometricControlPoints() const { return m_PhotometricControlPoints; }
	/// Sets the photometric control points.
	void SetPhotometricControlPoints(MatrixType const& photoCP) { m_PhotometricControlPoints = photoCP; }

	/// Gets the photometric weights.
	VectorType GetPhotometricWeights() const { return m_PhotometricWeights; }
	/// Sets the photometric weights.
	void SetPhotometricWeights(VectorType const& photoWeights) { m_PhotometricWeights = photoWeights; }


    ///	Return the type of the photometric kernel.
    inline KernelEnumType GetPhotometricKernelType() const { return m_PhotometricKernelType; }
    /// Set the type of the photometric kernel to \e kernelType.
    inline void SetPhotometricKernelType(KernelEnumType pkt) { m_PhotometricKernelType = pkt; }

    /// Gets the photometric gaussian kernel width.
    ScalarType GetPhotometricKernelWidth() const { return m_PhotometricKernelWidth; }
    /// Sets the photometric gaussian kernel width.
    void SetPhotometricKernelWidth(ScalarType const& pkw) { m_PhotometricKernelWidth = pkw; }

    /// Gets the photometric gaussian kernel field.
    MatrixType GetPhotometricKernelField() const { return m_PhotometricKernelField; }

    /// Gets the spacing parameter between photometric control points.
    ScalarType GetPhotometricCPSpacing() const { return m_PhotometricCPSpacing; }
    /// Sets the spacing parameter between photometric control points.
    void SetPhotometricCPSpacing(ScalarType const& spacing) { m_PhotometricCPSpacing = spacing; }


    /// Gets the number of photometric control points.
    unsigned long GetNumberOfPhotometricControlPoints() const { return m_NbPhotoCPs; }

    /// Returns the number of voxels in the down-sampled image domain.
    unsigned long GetNumberOfPoints() const {
        return m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels(); }
    /// Returns the dimension of the discretized image, here the number of voxels of the original image.
    unsigned long GetDimensionOfDiscretizedObject() const { return m_NbVoxels; }


    /// Returns the ITK reference image.
    inline ImagePointerType GetImage() const { return m_Image; }
    /// Returns the down-sampled ITK reference image.
    inline ImagePointerType GetDownSampledImage() const { return m_DownSampledImage; }
    /// Returns the coordinates of voxels in the domain of the down-sampled image.
    inline MatrixType GetDownSampledImageMap() const { return GridFunctionsType::ImageToPoints(m_DownSampledImage); }


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Other public method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

	/// Updates several information concerning the parametric image.
	virtual void Update();
    /// Saves in \e filename the parametric image.
    virtual void WriteObject(std::string filename) const;


	/// Returns the norm between itself and the target.
	virtual ScalarType ComputeMatch(const std::shared_ptr<Superclass> target);
	/// Returns the gradient of the norm between itself and the target.
	virtual MatrixType ComputeMatchGradient(const std::shared_ptr<Superclass> target);

    /// Returns the gradient of the data term wrt the photometric weights.
    MatrixType SplatDifferenceImage(const std::shared_ptr<LIImageType> target, const MatrixType& downSampledImageMap);


//    /// Computes the difference between the target and itself Ii - I, assuming a pixel-to-pixel correspondence.
//    VectorType Compute


    /// Initializes the photometric control points & weights, as well as itk-type "images".
    void Initialize(std::shared_ptr<LIImageType> const LIImage);
    /// Updates the photometric weights and related attributes.
    void UpdateImageIntensity(VectorType const& pw);


protected:

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Protected method(s) :
	////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Returns a deformed version of the image, where the deformation is given by the position of voxels positions in \e imagePoints.
    virtual std::shared_ptr<Superclass> doDeformedObject(const MatrixType& imagePoints) const {
        return std::static_pointer_cast<Superclass>(std::make_shared<ParametricImage>(*this, imagePoints)); }

    /// Makes a copy of the object.
    virtual std::shared_ptr<Superclass> doClone() const {
        return std::static_pointer_cast<Superclass>(std::make_shared<ParametricImage>(*this)); }


    /// Computes the bounding box of the image.
    void UpdateBoundingBox();

    /// Initializes the photometric control points as a regular lattice of spacing = m_PhotometricCPSpacing.
    void InitializePhotometricControlPoints();
    /// Initializes the photometric weights based on m_Image intensities at the control points.
    void InitializePhotometricWeights();
    /// Initializes the kernel field (convolution between m_Image voxel and m_PhotometricControlPoints positions).
    void InitializePhotometricKernelField(MatrixType const& grid);


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Attribute(s)
	////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Photometric control points ; size nbControlPoints x dimension.
    MatrixType m_PhotometricControlPoints;
    /// Photometric weights ; size nbControlPoints x 1.
    VectorType m_PhotometricWeights;

    /// Type of the photometric kernel.
    KernelEnumType m_PhotometricKernelType;
    /// Photometric kernel width.
    ScalarType m_PhotometricKernelWidth;
    /// Kernel scalar field on the image grid (m_NbPhotoCPs rows x m_NbVoxels cols).
    MatrixType m_PhotometricKernelField;

    /// Spacing parameter for photometric control points initialization.
    ScalarType m_PhotometricCPSpacing;
    /// Number of photometric control points.
    unsigned long m_NbPhotoCPs;
    /// Number of pixels / voxels in the reference image.
    unsigned long m_NbVoxels;
    
    /// Pointer to the reference ITK image and a down-sampled version.
    ImagePointerType m_Image;
    /// Pointer to a down-sampled version of the reference ITK image. Only voxel positions of this image are used.
    ImagePointerType m_DownSampledImage;
    /// Spatial gradient of the image, sampled on the m_Image grid.
    ImagePointerListType m_ImageSpatialGradient;
    
    /// Pointer to the reference LLImage, used for object writing.
    std::shared_ptr<LIImageType> m_LIImage;

    /// Vector indicating permutation among axes.
    std::vector<unsigned int> m_PermutationAxes;
    /// Vector indicating which axes are flipped (vector of 1s or -1s).
    std::vector<int> m_FlipAxes;
};
