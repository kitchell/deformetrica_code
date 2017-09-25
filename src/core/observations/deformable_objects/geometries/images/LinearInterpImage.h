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
#include "itkImage.h"
#include "GridFunctions.h"

#include <cstring>
#include <iostream>
#include <sstream>

/**
 *  \brief      An image using linear interpolation of voxel values.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The LinearInterpImage class inherited from AbstractGeometry models continuous images
 *              (i.e. with infinite resolution) by linearly interpolating image intensities within a pixel/voxel.\n\n
 *              A deformed image is obtained by setting voxels positions in space (\f$\phi^{-1}(v_k)\f$) and re-sampling
 *              image intensities at the original regular lattice of voxels (\f$v_k\f$).
 */
template<class ScalarType, unsigned int Dimension>
class LinearInterpImage : public AbstractGeometry<ScalarType, Dimension>
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
  /// List of ITK images
  typedef std::vector<ImageTypePointer> ImageList;

  /// Methods for computations with grids.
  typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  LinearInterpImage();
  virtual ~LinearInterpImage();

  /// Copy constructor.
  LinearInterpImage(const LinearInterpImage& other);
  /// Constructor which copies the object and resample the image.
  LinearInterpImage(const LinearInterpImage& example, const MatrixType& ImagePoints);

  /// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints.
  std::shared_ptr<LinearInterpImage> DeformedObject(const MatrixType& ImagePoints) const {
    return std::static_pointer_cast<LinearInterpImage>(doDeformedObject(ImagePoints)); }

  /// Clones the object.
  std::shared_ptr<LinearInterpImage> Clone() const { return std::static_pointer_cast<LinearInterpImage>(doClone()); }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets (ITK) image to \e img.
  void SetImage(ImageType* img) { SetImageAndDownSamplingFactor(img, 1); }
  // { m_Image = img; m_DownSampledY1.set_size(0,0);	this->SetModified(); }
  /// Sets (ITK) image and downsampling factor.
  virtual void SetImageAndDownSamplingFactor(ImageType* img, int ds);
  /// Updates intensities of the (ITK) image with \e I.
  void UpdateImageIntensity(const MatrixType& I);

  /// Returns the (ITK) image.
  ImageTypePointer GetImage() const { return m_Image; }
  /// Returns the (ITK) down-sampled image.
  ImageTypePointer GetDownSampledImage() const { return m_DownSampledImage; }
  /// Returns the list of image gradients.
  ImageList GetGradientImages() const { return m_GradientImages; }

  /// Returns the coordinates of voxels in the domain of the down-sampled image.
  MatrixType GetDownSampledImageMap() const { return GridFunctionsType::ImageToPoints(m_DownSampledImage); }
  // /// Set the coordinates of voxels in the domain of the down-sampled image.
  // void SetImagePoints(const MatrixType& Y1);


  /// Gets the permutation among axes.
  std::vector<unsigned int> GetPermutationAxes() const { return m_PermutationAxes; }
  /// Sets permutation among axes.
  void SetPermutationAxes(std::vector<unsigned int> permutation) { m_PermutationAxes = permutation; this->SetModified(); }

  /// Gets which axes have been flipped.
  std::vector<int> GetFlipAxes() { return m_FlipAxes; }
  /// Sets which axes have been flipped.
  void SetFlipAxes(std::vector<int> fa) { m_FlipAxes = fa; this->SetModified(); }


  /// Sets minimum of intensity range of the image (before rescaling).
  void SetMinIntensityOutput(double minI) { m_MinIntensityOutput = minI; }
  /// Sets maximum of intensity range of the image (before rescaling).
  void SetMaxIntensityOutput(double maxI) { m_MaxIntensityOutput = maxI; }

  // /// Return the number of voxels in the image domain (original, not down-sampled).
  // int GetNumberOfPoints() const { return m_NumberOfVoxels; }

  /// Returns the number of voxels in the down-sampled image domain.
  virtual unsigned long GetNumberOfPoints() const { return m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels(); }
  virtual unsigned long GetDimensionOfDiscretizedObject() const = 0;



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Update() = 0;

  virtual ScalarType ComputeMatch(const std::shared_ptr<Superclass> target) = 0;
  virtual MatrixType ComputeMatchGradient(const std::shared_ptr<Superclass> target) = 0;

  /// Splats residual image (\f$I_{source} - I_{target}\f$). Dual operation of interpolation.
  /// Used in particular for updating template image intensities in atlas construction.
  MatrixType SplatDifferenceImage(const std::shared_ptr<LinearInterpImage<ScalarType, Dimension>> target,
                                  const MatrixType& DownSampledImageMap);

  virtual void WriteObject(std::string filename) const;



 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns a Deformed version of the image, where the deformation is given by the position of voxels positions in \e ImagePoints.
  virtual std::shared_ptr<Superclass> doDeformedObject(const MatrixType& ImagePoints) const = 0;

  /// Clones the object.
  virtual std::shared_ptr<Superclass> doClone() const = 0;


  /// Computes the bounding box of the image.
  void UpdateBoundingBox();

  /// Upsample the mapping of the image points \e Y.
  MatrixType UpSampleImageMap(const MatrixType& Y);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Pointer to the reference ITK image.
  ImageTypePointer m_Image;
  /// Pointer to a down-sampled version of the reference ITK image. Only voxel positions of this image are used, not its intensities.
  ImageTypePointer m_DownSampledImage;
  /// List of gradient of ITK images.
  ImageList m_GradientImages;

  /// Minimum of image intensity range (before rescaling).
  double m_MinIntensityOutput;
  /// Maximum of image intensity range (before rescaling).
  double m_MaxIntensityOutput;
  /// Number of pixels/voxels in the reference image.
  unsigned long m_NumberOfVoxels;

  /// Vector indicating permutation among axes.
  std::vector<unsigned int> m_PermutationAxes;
  /// Vector indicating which axes are flipped (vector of 1s or -1s).
  std::vector<int> m_FlipAxes;


};

