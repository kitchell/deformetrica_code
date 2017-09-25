/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _GridFunctions_h
#define _GridFunctions_h

/// Support files.
#include "LinearAlgebra.h"

/// Librairies files.
#include <vector>
#include "itkImage.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkResampleImageFilter.h"
#include "LinearAlgebra.h"

using namespace def::algebra;

/**
 *  \brief      TODO .
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The GridFunctions class TODO .
 */
template<class ScalarType, unsigned int Dimension>
class GridFunctions {

 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// ITK image type.
  typedef itk::Image<ScalarType, Dimension> ImageType;
  /// ITK image pointer type.
  typedef typename ImageType::Pointer ImagePointer;
  /// ITK image index type.
  typedef typename ImageType::IndexType ImageIndexType;
  /// ITK image point type.
  typedef typename ImageType::PointType ImagePointType;
  /// TODO
  typedef typename ImageType::RegionType ImageRegionType;
  /// TODO
  typedef typename ImageType::SizeType ImageSizeType;
  /// TODO
  typedef typename ImageType::SpacingType ImageSpacingType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Convert image to set of values, ordered using standard ITK traversal order
  static VectorType VectorizeImage(const ImageType *img);

  /// Generate matrix of point locations from image grid.
  static MatrixType ImageToPoints(const ImageType *img);

  /// Returns a pointer on a (itk) image based on \e img whose the intensity of the voxels are given by \e values.
  static ImagePointer VectorToImage(const ImageType *img, const VectorType &values);

  /// Interpolate values from image at location \e pos.
  static VectorType Interpolate(const MatrixType &pos, const ImageType *values);

  /// TODO
  static ImagePointer SplatToImage(const ImageType *example, const MatrixType &pos,
                                   const VectorType &values, long gridPadding = 2);

  /// Downsample image to a size roughly old size / factor
  static ImagePointer DownsampleImage(const ImageType *img, ScalarType factor);

  /// Upsample \e img to the space defined by \e exImg.
  static ImagePointer UpsampleImage(const ImageType *exImg, const ImageType *img);

  /// TODO
  static MatrixType UpsampleImagePoints(const ImageType *img, const ImageType *downSampledImg,
                                        const MatrixType &downSampledPos);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Virtual method to avoid instanciation of the class GridFunctions.
  virtual void Abstract() = 0;

  /// TODO
  static inline void _getInterpolationWeightsAndGridPoints(std::vector<ScalarType> &weights,
                                                           std::vector<ImageIndexType> &gridIndices,
                                                           const VectorType &x,
                                                           const ImageType *img);

  /// Specialized linear interpolation weights for 2D.
  static inline void _getInterpolationWeightsAndGridPoints_2(std::vector<ScalarType> &weights,
                                                             std::vector<ImageIndexType> &gridIndices,
                                                             const VectorType &x,
                                                             const ImageType *img);

  /// Specialized linear interpolation weights for 3D.
  static inline void _getInterpolationWeightsAndGridPoints_3(std::vector<ScalarType> &weights,
                                                             std::vector<ImageIndexType> &gridIndices,
                                                             const VectorType &x,
                                                             const ImageType *img);

}; /* class GridFunctions */


#endif /* _GridFunctions_h */
