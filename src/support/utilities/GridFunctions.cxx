/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "GridFunctions.h"
#include <chrono>


////////////////////////////////////////////////////////////////////////////////////////////////////
// Public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
VectorType
GridFunctions<ScalarType, Dimension>
::VectorizeImage(const ImageType *img) {
  VectorType values(img->GetLargestPossibleRegion().GetNumberOfPixels(), 0);

  long r = 0;

  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  IteratorType it(img, img->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    values[r] = it.Get();

    r++;
  }

  return values;
}

template<class ScalarType, unsigned int Dimension>
MatrixType
GridFunctions<ScalarType, Dimension>
::ImageToPoints(const ImageType *img) {
  MatrixType pointsM(img->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension, 0);

  long r = 0;

  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  IteratorType it(img, img->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    ImageIndexType ind = it.GetIndex();

    ImagePointType p;
    img->TransformIndexToPhysicalPoint(ind, p);

    for (long dim = 0; dim < Dimension; dim++)
      pointsM(r, dim) = p[dim];
    r++;
  }

  return pointsM;
}

template<class ScalarType, unsigned int Dimension>
typename GridFunctions<ScalarType, Dimension>::ImagePointer
GridFunctions<ScalarType, Dimension>
::VectorToImage(const ImageType *imgEx, const VectorType &values) {
  if (imgEx->GetLargestPossibleRegion().GetNumberOfPixels() != values.size())
    throw std::runtime_error("Cannot set image voxels values: vector dimension mismatch");

  ImagePointer img = ImageType::New();
  img->SetRegions(imgEx->GetLargestPossibleRegion());
  img->CopyInformation(imgEx);
  img->Allocate();
  img->FillBuffer(0);

  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  IteratorType it(img, img->GetLargestPossibleRegion());

  unsigned int r = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    it.Set(values[r]);

    r++;
  }

  return img;
}

template<class ScalarType, unsigned int Dimension>
VectorType
GridFunctions<ScalarType, Dimension>
::Interpolate(const MatrixType &X, const ImageType *img) {
  long numPoints = X.rows();

  VectorType values(numPoints, 0.0);

  std::vector<ScalarType> weights;
  std::vector<ImageIndexType> gridIndices;

  for (unsigned int i = 0; i < numPoints; i++) {
    _getInterpolationWeightsAndGridPoints(weights, gridIndices, X.get_row(i), img);

    for (unsigned int j = 0; j < weights.size(); j++) {
      ScalarType w = weights[j];
      ImageIndexType ind = gridIndices[j];
      values[i] += w * img->GetPixel(ind);
    }
  }

  return values;
}

template<class ScalarType, unsigned int Dimension>
typename GridFunctions<ScalarType, Dimension>::ImagePointer
GridFunctions<ScalarType, Dimension>
::SplatToImage(const ImageType *exImg, const MatrixType &X, const VectorType &values, long gridPadding) {
  long numPoints = X.rows();

  ImagePointer img = ImageType::New();
  img->SetRegions(exImg->GetLargestPossibleRegion());
  img->CopyInformation(exImg);
  img->Allocate();
  img->FillBuffer(0);

  ImageSizeType size = exImg->GetLargestPossibleRegion().GetSize();

  for (unsigned int i = 0; i < numPoints; i++) {
    std::vector<ScalarType> weights;
    std::vector<ImageIndexType> gridIndices;

    _getInterpolationWeightsAndGridPoints(weights, gridIndices, X.get_row(i), img);

    for (unsigned int j = 0; j < weights.size(); j++) {
      ScalarType w = weights[j];
      ImageIndexType ind = gridIndices[j];

      bool isborder = false;
      for (unsigned int d = 0; d < Dimension; d++)
        if (ind[d] <= gridPadding || ind[d] >= ((long) size[d] - gridPadding))
          isborder = true;
      if (isborder)
        continue;

      img->SetPixel(ind, img->GetPixel(ind) + w * values[i]);
    }
  }
  return img;
}

template<class ScalarType, unsigned int Dimension>
typename GridFunctions<ScalarType, Dimension>::ImagePointer
GridFunctions<ScalarType, Dimension>
::DownsampleImage(const ImageType *img, ScalarType factor) {
  if (factor <= 1.0) {
    typedef itk::ImageDuplicator<ImageType> DuperType;
    typename DuperType::Pointer dupef = DuperType::New();
    dupef->SetInputImage(img);
    dupef->Update();
    return dupef->GetOutput();
  }

  ImageSpacingType spacing = img->GetSpacing();
  ImageSizeType size = img->GetLargestPossibleRegion().GetSize();

  ScalarType minSpacing = spacing[0];
  for (unsigned int dim = 1; dim < Dimension; dim++)
    if (spacing[dim] < minSpacing)
      minSpacing = spacing[dim];

  ImageSpacingType downSpacing;
  ImageSizeType downSize;

  // MARCEL
  // for (unsigned int dim = 0; dim < Dimension; dim++)
  // {
  //   ScalarType length = size[dim]*spacing[dim];
  //
  //   //downSize[dim] = static_cast<long>(length / factor);
  //   downSize[dim] = static_cast<long>(size[dim] / factor);
  //
  //   downSpacing[dim] = length / downSize[dim];
  // }

  // STANLEY: this function is used to define a coarser voxel grid, and then to interpolate it to give the coordinate of voxels at a finest resolution
  // in this case, the downsampled image needs to cover a *larger* domain than the fine image. Otherwise, voxels in the fine grid outside the domain covered by the coarse grid will be interpolated with zero values in the neighborhood located outside the domain
  // If this is desirable to interpolate intensities (assuming black background), this has dramatic consequences when interpolating voxel coordinates


  for (unsigned int dim = 0; dim < Dimension; dim++) {
    downSize[dim] = ceil(size[dim] / factor) + 1;
    downSpacing[dim] = spacing[dim] * factor;
  }

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>
      GaussFilterType;
  typename GaussFilterType::Pointer gaussf = GaussFilterType::New();
  gaussf->SetInput(img);
  gaussf->SetVariance(minSpacing * minSpacing);
  //gaussf->SetMaximumKernelWidth(5);
  //gaussf->SetMaximumError(0.1);
  gaussf->Update();

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
  typename ResamplerType::Pointer resf = ResamplerType::New();
  //resf->SetInput(img);
  resf->SetInput(gaussf->GetOutput());
  resf->SetOutputParametersFromImage(img);
  resf->SetSize(downSize);
  resf->SetOutputSpacing(downSpacing);
  resf->Update();

  return resf->GetOutput();
}

template<class ScalarType, unsigned int Dimension>
typename GridFunctions<ScalarType, Dimension>::ImagePointer
GridFunctions<ScalarType, Dimension>
::UpsampleImage(const ImageType *exImg, const ImageType *img) {
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
  typename ResamplerType::Pointer resf = ResamplerType::New();
  resf->SetInput(img);
  resf->SetOutputParametersFromImage(exImg);
  resf->Update();

  return resf->GetOutput();
}

template<class ScalarType, unsigned int Dimension>
MatrixType
GridFunctions<ScalarType, Dimension>
::UpsampleImagePoints(const ImageType *img, const ImageType *downSampledImg, const MatrixType &downSampledPos) {
  unsigned int numPixels = img->GetLargestPossibleRegion().GetNumberOfPixels();

  MatrixType Yup(numPixels, Dimension, 0.0);

  MatrixType PointImage = GridFunctions<ScalarType, Dimension>::ImageToPoints(img);
  for (unsigned int d = 0; d < Dimension; d++) {
    ImagePointer Hd = GridFunctions<ScalarType, Dimension>::VectorToImage(downSampledImg, downSampledPos.get_column(d));
    Yup.set_column(d, GridFunctions<ScalarType, Dimension>::Interpolate(PointImage, Hd));
  }

  return Yup;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class ScalarType, unsigned int Dimension>
void
GridFunctions<ScalarType, Dimension>
::_getInterpolationWeightsAndGridPoints(std::vector<ScalarType> &weights,
                                        std::vector<ImageIndexType> &gridIndices,
                                        const VectorType &x,
                                        const ImageType *img) {
  if (Dimension == 2) {
    // Bilinear
    _getInterpolationWeightsAndGridPoints_2(weights, gridIndices, x, img);
  } else if (Dimension == 3) {
    // Trilinear
    _getInterpolationWeightsAndGridPoints_3(weights, gridIndices, x, img);
  } else {
    // Nearest neighbor
    ImagePointType p;
    for (unsigned int k = 0; k < Dimension; k++)
      p[k] = x[k];

    ImageSizeType size = img->GetLargestPossibleRegion().GetSize();

    // we need to make sure that this method well behaves when point coordinates are negative (floor instead of int conversion)
    ImageIndexType ind;
    img->TransformPhysicalPointToIndex(p, ind);

    weights.clear();
    gridIndices.clear();

    bool isOut = false;
    for (unsigned int k = 0; k < Dimension; k++)
      if (ind[k] < 0 || ind[k] >= (long) size[k]) {
        isOut = true;
        break;
      }

    if (!isOut) {
      weights.push_back(1.0);
      gridIndices.push_back(ind);
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
GridFunctions<ScalarType, Dimension>
::_getInterpolationWeightsAndGridPoints_2(std::vector<ScalarType> &weights,
                                          std::vector<ImageIndexType> &gridIndices,
                                          const VectorType &p,
                                          const ImageType *img) {
  //  ImagePointType origin = img->GetOrigin();
  //  ImageSpacingType spacing = img->GetSpacing();

  ImageSizeType size = img->GetLargestPossibleRegion().GetSize();

  // MARCEL
  // std::vector<ScalarType> contInd(Dimension);
  // for (unsigned int d = 0; d < Dimension; d++)
  //   contInd[d] = (p[d]-origin[d]) / spacing[d];

  //STANLEY: solve problem when axes are flipped and transform contInd[d] = (p[d]-origin[d]) / spacing[d] is inaccurate
  typedef itk::Point<ScalarType, Dimension> PointType;
  typedef itk::ContinuousIndex<ScalarType, Dimension> ContIndexType;
  PointType pt;
  for (unsigned int d = 0; d < Dimension; d++)
    pt[d] = p[d];

  ContIndexType contInd;
  img->TransformPhysicalPointToContinuousIndex(pt, contInd);


  //
  // Bilinear interpolation
  //

  // Get the 4 grid positions
  // STANLEY: floor and not int conversion, since contIndex may be negative (in case deformation field ends up outside image boundary for instance)
  // but since it may also have a neighborhing voxel *on* the boundary, one wants to keep it
  int ix1 = floor(contInd[0]);
  int iy1 = floor(contInd[1]);

  int ix2 = ix1 + 1;
  int iy2 = iy1 + 1;

  ScalarType fx = contInd[0] - ix1;
  ScalarType fy = contInd[1] - iy1;

  ScalarType gx = ix2 - contInd[0];
  ScalarType gy = iy2 - contInd[1];

  // Add valid grid positions and corresponding weights
  weights.clear();
  gridIndices.clear();

#define interpWeightMacro2(ix, iy, w) \
        if ((0 <= (ix)) && ((ix) < (int)size[0]) && \
                (0 <= (iy)) && ((iy) < (int)size[1])) \
                { \
            ImageIndexType ind; \
            ind[0] = (ix); ind[1] = (iy); \
            if (w > 0) \
            { \
                gridIndices.push_back(ind); \
                weights.push_back(w); \
            } \
                }

  interpWeightMacro2(ix1, iy1, gx * gy);
  interpWeightMacro2(ix1, iy2, gx * fy);
  interpWeightMacro2(ix2, iy1, fx * gy);
  interpWeightMacro2(ix2, iy2, fx * fy);

#undef interpWeightMacro2
}

template<class ScalarType, unsigned int Dimension>
void
GridFunctions<ScalarType, Dimension>
::_getInterpolationWeightsAndGridPoints_3(std::vector<ScalarType> &weights,
                                          std::vector<ImageIndexType> &gridIndices,
                                          const VectorType &p,
                                          const ImageType *img) {
  //  ImagePointType origin = img->GetOrigin();
  //  ImageSpacingType spacing = img->GetSpacing();

  ImageSizeType size = img->GetLargestPossibleRegion().GetSize();


  // MARCEL
  //   std::vector<ScalarType> contInd(Dimension);
  //   for (unsigned int d = 0; d < Dimension; d++)
  //   {
  // 	contInd[d] = (p[d]-origin[d]) / spacing[d];
  // }

  //STANLEY: solve problem when axes are flipped and transform contInd[d] = (p[d]-origin[d]) / spacing[d] is inaccurate
  typedef itk::Point<ScalarType, Dimension> PointType;
  typedef itk::ContinuousIndex<ScalarType, Dimension> ContIndexType;
  PointType pt;
  for (unsigned int d = 0; d < Dimension; d++)
    pt[d] = p[d];

  ContIndexType contInd;
  img->TransformPhysicalPointToContinuousIndex(pt, contInd);

  //
  // Trilinear interpolation
  //

  // Get the 8 grid positions
  // STANLEY: floor and not int conversion, since contIndex may be negative (in case deformation field ends up outside image boundary for instance)
  // but since it may also have a neighborhing voxel *on* the boundary, one wants to keep it
  int ix1 = floor(contInd[0]);
  int iy1 = floor(contInd[1]);
  int iz1 = floor(contInd[2]);

  int ix2 = ix1 + 1;
  int iy2 = iy1 + 1;
  int iz2 = iz1 + 1;

  // Get distances to the image grid
  ScalarType fx = contInd[0] - ix1;
  ScalarType fy = contInd[1] - iy1;
  ScalarType fz = contInd[2] - iz1;

  ScalarType gx = ix2 - contInd[0];
  ScalarType gy = iy2 - contInd[1];
  ScalarType gz = iz2 - contInd[2];

  // Add valid grid positions and corresponding weights
  weights.clear();
  gridIndices.clear();

#define interpWeightMacro3(ix, iy, iz, w) \
        if ((0 <= (ix)) && ((ix) < (int)size[0]) && \
                (0 <= (iy)) && ((iy) < (int)size[1]) && \
                (0 <= (iz)) && ((iz) < (int)size[2])) \
                { \
            ImageIndexType ind; \
            ind[0] = (ix); ind[1] = (iy); ind[2] = (iz); \
            if (w > 0) \
            { \
                gridIndices.push_back(ind); \
                weights.push_back(w); \
            } \
                }

  interpWeightMacro3(ix1, iy1, iz1, gx * gy * gz);
  interpWeightMacro3(ix1, iy1, iz2, gx * gy * fz);
  interpWeightMacro3(ix1, iy2, iz1, gx * fy * gz);
  interpWeightMacro3(ix1, iy2, iz2, gx * fy * fz);
  interpWeightMacro3(ix2, iy1, iz1, fx * gy * gz);
  interpWeightMacro3(ix2, iy1, iz2, fx * gy * fz);
  interpWeightMacro3(ix2, iy2, iz1, fx * fy * gz);
  interpWeightMacro3(ix2, iy2, iz2, fx * fy * fz);

#undef interpWeightMacro3
}

template class GridFunctions<double, 2>;
template class GridFunctions<double, 3>;
