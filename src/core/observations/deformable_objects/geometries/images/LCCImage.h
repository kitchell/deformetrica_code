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

#include "LinearInterpImage.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      An image using linear interpolation of voxel values and LCC Metric between images.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The LCCImage class inherited from LinearInterpImage ... TODO .
 */
template<class ScalarType, unsigned int Dimension>
class LCCImage : public LinearInterpImage<ScalarType, Dimension>
{
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract image type.
  typedef LinearInterpImage<ScalarType, Dimension> Superclass;

  /// Deformable Object type.
  typedef typename Superclass::Superclass AbstractGeometryType;

  /// ITK image type.
  typedef itk::Image<ScalarType, Dimension> ImageType;
  /// ITK image pointer type.
  typedef typename ImageType::Pointer ImageTypePointer;
  /// ITK image index type.
  typedef typename ImageType::IndexType ImageIndexType;
  /// ITK image point type.
  typedef typename ImageType::PointType ImagePointType;
  /// ITK image region type.
  typedef typename ImageType::RegionType ImageRegionType;
  /// ITK image size type.
  typedef typename ImageType::SizeType ImageSizeType;
  /// ITK image spacing type.
  typedef typename ImageType::SpacingType ImageSpacingType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  LCCImage();

  /// Destructor.
  virtual ~LCCImage();

  /// Copy constructor.
  LCCImage(const LCCImage& other);
  /// Constructor which copies the object and resample the image
  LCCImage(const LCCImage& example, const MatrixType& ImagePoints);

  /// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints
  std::shared_ptr<LCCImage> DeformedObject(const MatrixType& ImagePoints) const {
	  return std::static_pointer_cast<LCCImage>(doDeformedObject(ImagePoints)); }

  /// Clones the object.
  std::shared_ptr<LCCImage> Clone() const {
	  return std::static_pointer_cast<LCCImage>(doClone()); }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Set standard deviation of the Gaussian filter needed to compute local statistics
  void SetLCCKernelWidth(ScalarType h) { m_LCCKernelWidth = h; }

  /// Get Local Mean Image
  ImageTypePointer GetLocalMeanImage() const {return m_LocalMeanImage; }
  /// Get Local Standard Deviation Image
  ImageTypePointer GetLocalVarianceImage() const { return m_LocalVarianceImage; }

  virtual unsigned long GetDimensionOfDiscretizedObject() const { throw std::runtime_error("No noise model for LCC Image!"); return 0; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  void Update();

  /// Compute LCC between this image (after deformation) and target image
  virtual ScalarType ComputeMatch(const std::shared_ptr<AbstractGeometryType> target);

  /// Compute the gradient of the LCC metric between this image (once deformed) and target image
  virtual MatrixType ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target);


 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints
  virtual std::shared_ptr<AbstractGeometryType> doDeformedObject(const MatrixType& ImagePoints) const {
	  return std::static_pointer_cast<AbstractGeometryType>(std::make_shared<LCCImage>(*this, ImagePoints)); }

  /// Clones the object.
  virtual std::shared_ptr<AbstractGeometryType> doClone() const {
	  return std::static_pointer_cast<AbstractGeometryType>(std::make_shared<LCCImage>(*this)); }


  /// Computes the local mean image from input image, namely the convolution between input image and a gaussian filter: W*I
  ImageTypePointer ComputeLocalMeanImage(const ImageType* img) const;

  /// Computes the local covariance image from 2 input images (requires to have Local Mean Image computed): W*(I.J) - (W*I)(W*J)
  ImageTypePointer ComputeLocalCovarianceImage(const ImageType* img1, const ImageType* LocalMeanImg1,
											   const ImageType* img2, const ImageType* LocalMeanImg2, bool sameImage) const;



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Image of local mean (convolution of the image I with kernel W: W*I)
  ImageTypePointer m_LocalMeanImage;

  /// Image of local standard deviation: sqrt( (W*(I^2)) )
  ImageTypePointer m_LocalVarianceImage;


  /// Standard deviation of the Gaussian function used to compute the local mean and correlation
  ScalarType m_LCCKernelWidth;

  /// Minimum image size
  // int m_SmallestImageDimension;


};
