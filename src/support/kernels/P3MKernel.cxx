/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/


#include "P3MKernel.h"

#include "itkConfigure.h"
#if defined(USE_FFTWF)
#if ITK_VERSION_MAJOR < 4
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"
#else
#include "itkFFTWForwardFFTImageFilter.h"
#include "itkFFTWInverseFFTImageFilter.h"
#endif
#else
#if ITK_VERSION_MAJOR < 4
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#else
#include "itkVnlForwardFFTImageFilter.h"
#include "itkVnlInverseFFTImageFilter.h"
#endif
#endif

#include "itkVersion.h"

#include <itkFFTShiftImageFilter.h>

#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"

#include <exception>
#include <stdexcept>

#define DO_NEAR_FIELD 0
#define DO_SORT_ID 0

//
// Support class
//

class TagIndex {
 public:
  unsigned int index;
  unsigned int j;
  inline TagIndex &operator=(const TagIndex &o) {
    this->index = o.index;
    this->j = o.j;
    return (*this);
  }
  inline bool operator<(const TagIndex &o) const { return this->index < o.index; }
};

//
// P3MKernel class definition
//

// Static variables
template<class ScalarType, unsigned int PointDim>
itk::SimpleFastMutexLock
    P3MKernel<ScalarType, PointDim>::m_FFTMutex;

template<class ScalarType, unsigned int PointDim>
itk::SimpleFastMutexLock
    P3MKernel<ScalarType, PointDim>::m_CacheMutex;

template<class ScalarType, unsigned int PointDim>
std::vector<ScalarType>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTKernelWidths;

template<class ScalarType, unsigned int PointDim>
std::vector<typename P3MKernel<ScalarType, PointDim>::ImageSizeType>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTKernelSizes;

template<class ScalarType, unsigned int PointDim>
std::vector<typename P3MKernel<ScalarType, PointDim>::ComplexImagePointer>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTKernelImages;

template<class ScalarType, unsigned int PointDim>
std::vector<ScalarType>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTGradientKernelWidths;

template<class ScalarType, unsigned int PointDim>
std::vector<typename P3MKernel<ScalarType, PointDim>::ImageSizeType>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTGradientKernelSizes;

template<class ScalarType, unsigned int PointDim>
std::vector<std::vector<typename P3MKernel<ScalarType, PointDim>::ComplexImagePointer> >
    P3MKernel<ScalarType, PointDim>::m_CacheFFTGradientKernelImages;

template<class ScalarType, unsigned int PointDim>
std::vector<ScalarType>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTHessianKernelWidths;

template<class ScalarType, unsigned int PointDim>
std::vector<typename P3MKernel<ScalarType, PointDim>::ImageSizeType>
    P3MKernel<ScalarType, PointDim>::m_CacheFFTHessianKernelSizes;

template<class ScalarType, unsigned int PointDim>
std::vector<std::vector<typename P3MKernel<ScalarType, PointDim>::ComplexImagePointer> >
    P3MKernel<ScalarType, PointDim>::m_CacheFFTHessianKernelImages;


///// For profiling.
//template<class ScalarType, unsigned int PointDim> ScalarType P3MKernel<ScalarType, PointDim>::m_ConvolveTime;
//template<class ScalarType, unsigned int PointDim> ScalarType P3MKernel<ScalarType, PointDim>::m_ConvolveGradientTime;
//template<class ScalarType, unsigned int PointDim> ScalarType P3MKernel<ScalarType, PointDim>::m_ConvolveHessianTime;
//
//template<class ScalarType, unsigned int PointDim> ScalarType P3MKernel<ScalarType, PointDim>::m_UpdateGridTime;
//template<class ScalarType, unsigned int PointDim> ScalarType P3MKernel<ScalarType, PointDim>::m_ApplyKernelFftTime;
//template<class ScalarType, unsigned int PointDim> ScalarType P3MKernel<ScalarType, PointDim>::m_InterpolateTime;


// Methods

template<class ScalarType, unsigned int PointDim>
P3MKernel<ScalarType, PointDim>
::P3MKernel(): Superclass() {
  this->Init();
  this->SetKernelWidth(1.0);
}

template<class ScalarType, unsigned int PointDim>
P3MKernel<ScalarType, PointDim>
::P3MKernel(const P3MKernel &o) {
  //this->SetSources(o.GetSources());
  //this->SetWeights(o.GetWeights());
  Superclass::m_Sources = o.m_Sources;
  Superclass::m_Weights = o.m_Weights;
  this->SetKernelWidth(o.GetKernelWidth());

  if (o.IsModified())
    this->SetModified();
  else
    this->UnsetModified();

  m_GridOrigin = o.m_GridOrigin;
  m_GridSpacing = o.m_GridSpacing;
  m_GridSize = o.m_GridSize;
  m_GridPadding = o.m_GridPadding;

  m_FFTKernel = o.m_FFTKernel;
  m_FFTGradientKernels = o.m_FFTGradientKernels;
  m_FFTHessianKernels = o.m_FFTHessianKernels;

  // m_MeshPreConvList  = o.m_MeshPreConvList;

  m_MeshList = o.m_MeshList;
  // m_MeshWYList  = o.m_MeshWYList;
  // m_MeshWYYtList  = o.m_MeshWYYtList;

  m_MeshGradientList = o.m_MeshGradientList;
  m_MeshHessianList = o.m_MeshHessianList;

  m_HessianUpdated = o.m_HessianUpdated;

  m_SourceIDs = o.m_SourceIDs;

  m_IDLookupTable = o.m_IDLookupTable;

  m_NearThresholdScale = o.m_NearThresholdScale;
}

template<class ScalarType, unsigned int PointDim>
P3MKernel<ScalarType, PointDim>
::P3MKernel(const MatrixType &X, double h): Superclass(X, h) {
  this->Init();
  this->SetKernelWidth(h);
  this->SetModified();
}

template<class ScalarType, unsigned int PointDim>
P3MKernel<ScalarType, PointDim>
::P3MKernel(
    const MatrixType &X, const MatrixType &W,
    double h): Superclass(X, W, h) {
  this->Init();
  this->SetKernelWidth(h);
  this->SetModified();
}

template<class ScalarType, unsigned int PointDim>
P3MKernel<ScalarType, PointDim>
::~P3MKernel() {
  this->ClearGrids();
  m_FFTKernel = 0;
}

template<class ScalarType, unsigned int PointDim>
P3MKernel<ScalarType, PointDim> *
P3MKernel<ScalarType, PointDim>
::Clone() const {
  return new P3MKernel(*this);
}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::ClearGrids() {

  m_MeshList.clear();
  m_MeshGradientList.clear();
  m_MeshHessianList.clear();
}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::Init() {
  m_GridOrigin.Fill(0.0);
  m_GridSpacing.Fill(1.0);
  m_GridSize.Fill(128);
  m_GridPadding = 0;
  m_NearThresholdScale = 0.5f;
  m_FFTKernel = 0;

  m_HessianUpdated = true;

  m_WorkingSpacingRatio = 1.0;
  m_PaddingFactor = 0.0;

}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::UpdateGrids() {

  // Update grid size and spacing
  this->DetermineGrids();

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  ImageRegionType region;
  region.SetSize(m_GridSize);

  this->ClearGrids();

  for (unsigned int k = 0; k < weightDim; k++) {
    ImagePointer img = ImageType::New();
    img->SetRegions(region);
    img->Allocate();
    img->SetOrigin(m_GridOrigin);
    img->SetSpacing(m_GridSpacing);
    img->FillBuffer(0);

    m_MeshList.push_back(img);
  }

  // Splat weight values to meshes
  for (unsigned int i = 0; i < Superclass::m_Sources.rows(); i++) {
    VectorType yi = Superclass::m_Sources.get_row(i);

    std::vector<ScalarType> weights;
    std::vector<ImageIndexType> gridIndices;

    this->_getInterpolationWeightsAndGridPoints(
        weights, gridIndices, yi);

    for (unsigned int j = 0; j < weights.size(); j++) {
      ScalarType w = weights[j];
      ImageIndexType ind = gridIndices[j];

      bool isborder = false;
      for (unsigned int d = 0; d < PointDim; d++)
        if (ind[d] <= m_GridPadding || ind[d] >= ((long) m_GridSize[d] - m_GridPadding))
          isborder = true;
      if (isborder)
        continue;

      for (unsigned int k = 0; k < weightDim; k++)
        m_MeshList[k]->SetPixel(ind,
                                m_MeshList[k]->GetPixel(ind) + w * Superclass::m_Weights(i, k));
    }
  }

  // apply FFT to the weight meshes
#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::FFTWComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::FFTWInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::VnlInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#endif

  // TODO:
  // fftw execute is thread safe, but planner is not
  // for thread safety and speed may need to write own ITK fftw wrapper

  m_MeshListFFT.resize(weightDim);
  for (unsigned int k = 0; k < weightDim; k++) {
    m_MeshList[k]->ReleaseDataFlagOff();

    //P3MKernel::m_FFTMutex.Lock();

    typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
    forwardFFT->SetInput(m_MeshList[k]);
    forwardFFT->SetNumberOfThreads(1);
    forwardFFT->Update();

    //P3MKernel::m_FFTMutex.Unlock();

    m_MeshListFFT[k] = forwardFFT->GetOutput();
  }

  m_FFTKernel = this->BuildFFTKernel();
  m_FFTGradientKernels = this->BuildFFTGradientKernels();

  m_HessianUpdated = false;
}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::UpdateHessianGrids() {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
//	unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Source and weights count mismatch");


  m_FFTHessianKernels = this->BuildFFTHessianKernels();

  m_HessianUpdated = true;

}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::DetermineGrids() {
  VectorType Xmin = m_DataDomain.get_column(0);
  VectorType Xmax = m_DataDomain.get_column(1);

  ScalarType h = this->GetKernelWidth();

  if (h <= 1e-20) {
    //     this->SetGridOrigin(origin);
    //     this->SetGridSize(size);
    //     this->SetGridSpacing(spacing);
    // this->SetGridPadding(0.0);
    std::cout << "kernel width too small to set up a p3m grid!" << std::endl;
    return;
  }

  ImageSpacingType grid_spacing;
  ImageSizeType grid_size;
  ImagePointType grid_origin;

  grid_spacing.Fill(m_WorkingSpacingRatio * h);

  ScalarType length;
  ScalarType padded_length;
  for (unsigned int d = 0; d < PointDim; d++) {
    length = Xmax[d] - Xmin[d]; //size[d] * spacing[d];
    padded_length = length + 2 * m_PaddingFactor * h;

    grid_size[d] = (long) (padded_length / grid_spacing[d]);
  }

  std::vector<long> pads(PointDim, 0);
  for (unsigned int d = 0; d < PointDim; d++) {
    // log_2(size)
    ScalarType l2 = log((ScalarType) grid_size[d]) / log(2.0);

    // Round up
    long trunc_l2 = (long) l2;
    if ((l2 - trunc_l2) > 0)
      l2 = trunc_l2 + 1;
    else
      l2 = trunc_l2;

    // Nearest power of 2
    grid_size[d] = (long) pow(2, l2);

    ScalarType shift = (grid_size[d] * grid_spacing[d] - (Xmax[d] - Xmin[d]))
        / 2; //(grid_size[d]*grid_spacing[d] - size[d]*spacing[d]) / 2;
    grid_origin[d] = Xmin[d] - shift; //origin[d] - shift;

  }

  // if (grid_size != m_GridSize)
  // std::cout << "Grid changed: origin = " << grid_origin << " size = " << grid_size << " spacing = " << grid_spacing << std::endl;


  this->SetGridOrigin(grid_origin);
  this->SetGridSize(grid_size);
  this->SetGridSpacing(grid_spacing);

}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::_getInterpolationWeightsAndGridPoints(
    std::vector<ScalarType> &weights,
    std::vector<ImageIndexType> &gridIndices,
    const VectorType &x) {
  //TODO sync with the one in GridFunctions
  //typedef GridFunctions<ScalarType, PointDim> GridFunctionsType;
  //GridFunctionsType::_getInterpolationWeightsAndGridPoints(weights, gridIndices, x, m_MeshList[0]);

  if (PointDim == 2) {
    // Bilinear
    _getInterpolationWeightsAndGridPoints_2(weights, gridIndices, x);
  } else if (PointDim == 3) {
    // Trilinear
    _getInterpolationWeightsAndGridPoints_3(weights, gridIndices, x);
  } else {
    // Nearest neighbor
    ImageIndexType ind;
    for (unsigned int d = 0; d < PointDim; d++)
      ind[d] = (long) ((x[d] - m_GridOrigin[d]) / m_GridSpacing[d]);
    gridIndices.clear();

    bool isOut = false;
    for (unsigned int k = 0; k < PointDim; k++)
      if (ind[k] < 0 || ind[k] >= (long) m_GridSize[k]) {
        isOut = true;
        break;
      }

    if (!isOut) {
      weights.push_back(1.0);
      gridIndices.push_back(ind);
    }
  }
}

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::_getInterpolationWeightsAndGridPoints_2(
    std::vector<ScalarType> &weights,
    std::vector<ImageIndexType> &gridIndices,
    const VectorType &p) {
  // it works because origin is always the point with the smallest coordinate in each direction (i.e. axes are not flipped) -- different from same function in GridFunctions.txx
  std::vector<ScalarType> contInd(PointDim);
  for (unsigned int d = 0; d < PointDim; d++)
    contInd[d] = (p[d] - m_GridOrigin[d]) / m_GridSpacing[d];

  //
  // Bilinear interpolation
  //

  // Get the 4 grid positions

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
  if ((0 <= (ix)) && ((ix) < (int)m_GridSize[0]) && \
    (0 <= (iy)) && ((iy) < (int)m_GridSize[1])) \
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

template<class ScalarType, unsigned int PointDim>
void
P3MKernel<ScalarType, PointDim>
::_getInterpolationWeightsAndGridPoints_3(
    std::vector<ScalarType> &weights,
    std::vector<ImageIndexType> &gridIndices,
    const VectorType &p) {
  // it works because origin is always the point with the smallest coordinate in each direction (i.e. axes are not flipped) -- different from same function in GridFunctions.txx
  std::vector<ScalarType> contInd(PointDim);
  for (unsigned int d = 0; d < PointDim; d++)
    contInd[d] = (p[d] - m_GridOrigin[d]) / m_GridSpacing[d];

  //
  // Trilinear interpolation
  //

  // Get the 8 grid positions
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
  if ((0 <= (ix)) && ((ix) < (int)m_GridSize[0]) && \
    (0 <= (iy)) && ((iy) < (int)m_GridSize[1]) && \
    (0 <= (iz)) && ((iz) < (int)m_GridSize[2])) \
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

template<class ScalarType, unsigned int PointDim>
long
P3MKernel<ScalarType, PointDim>
::_getPointID(const VectorType &x) {
  std::vector<ScalarType> contInd(PointDim);
  for (unsigned int d = 0; d < PointDim; d++)
    contInd[d] = (x[d] - m_GridOrigin[d]) / m_GridSpacing[d];

  // Check if out of bounds
  bool isOut = false;
  for (unsigned int d = 0; d < PointDim; d++)
    if (contInd[d] < 0 || contInd[d] > (m_GridSize[d] - 1)) {
      isOut = true;
      break;
    }

  if (isOut) {
    long outid = 1;
    for (unsigned int d = 0; d < PointDim; d++)
      outid *= m_GridSize[d];
    outid += 1;
    return outid;
  }

  // Round up continiuous index to nearest point
  ImageIndexType ind;
  for (unsigned int d = 0; d < PointDim; d++) {
    ind[d] = (long) (contInd[d] + 0.5f);
  }

  return this->_getPointID(ind);
}

template<class ScalarType, unsigned int PointDim>
long
P3MKernel<ScalarType, PointDim>
::_getPointID(const ImageIndexType &ind) {
  long id = ind[PointDim - 1];
  long mult = m_GridSize[PointDim - 1];
  for (long d = PointDim - 2; d >= 0; d--) {
    id += ind[d] * mult;
    mult *= m_GridSize[d];
  }

  return id;
}

template<class ScalarType, unsigned int PointDim>
VectorType
P3MKernel<ScalarType, PointDim>
::Interpolate(
    const VectorType &x, const std::vector<ImagePointer> &images) {
  unsigned int numInterps = images.size();

  VectorType values(numInterps, 0.0);

  std::vector<ScalarType> weights;
  std::vector<ImageIndexType> gridIndices;

  this->_getInterpolationWeightsAndGridPoints(
      weights, gridIndices, x);

  for (unsigned int i = 0; i < weights.size(); i++) {
    ScalarType w = weights[i];
    ImageIndexType ind = gridIndices[i];
    for (unsigned int j = 0; j < numInterps; j++)
      values[j] += w * images[j]->GetPixel(ind);
  }

  return values;
}

template<class ScalarType, unsigned int PointDim>
typename P3MKernel<ScalarType, PointDim>::ComplexImagePointer
P3MKernel<ScalarType, PointDim>
::BuildFFTKernel() {
  ScalarType kernelWidth = this->GetKernelWidth();

  m_CacheMutex.Lock();

  unsigned int numCache = m_CacheFFTKernelWidths.size();

  unsigned int cacheIndex = numCache + 1;
  for (unsigned int i = 0; i < numCache; i++)
    if (m_CacheFFTKernelWidths[i] == kernelWidth) {
      bool sameSize = true;
      ImageSizeType size_i = m_CacheFFTKernelSizes[i];
      for (unsigned int d = 0; d < PointDim; d++)
        if (size_i[d] != m_GridSize[d]) {
          sameSize = false;
          break;
        }
      if (sameSize) {
        cacheIndex = i;
        break;
      }
    }

  m_CacheMutex.Unlock();

  if (cacheIndex < numCache) {
    m_CacheMutex.Lock();
    ComplexImagePointer kImg = m_CacheFFTKernelImages[cacheIndex];
    m_CacheMutex.Unlock();

    return kImg;
  } else {
    std::cout << "New grid size: " << m_GridSize << " with spacing = " << m_GridSpacing << std::endl;
  }


//std::cout << "Building kernel with width = " << kernelWidth << " and size = " << m_GridSize << std::endl;

  ImagePointer kernelImg = ImageType::New();
  kernelImg->SetRegions(m_MeshList[0]->GetLargestPossibleRegion());
  kernelImg->Allocate();
  kernelImg->CopyInformation(m_MeshList[0]);
  kernelImg->FillBuffer(0);

  ScalarType nearThres = 0;
  for (unsigned int d = 0; d < PointDim; d++)
    nearThres += m_GridSpacing[d] * m_GridSpacing[d];
  nearThres *= m_NearThresholdScale * m_NearThresholdScale;

  VectorType center(PointDim, 0.0);

  for (unsigned int d = 0; d < PointDim; d++) {
    center[d] =
        m_GridOrigin[d] + 0.5f * m_GridSize[d] * m_GridSpacing[d];
  }

  ImageIteratorType it(
      kernelImg, kernelImg->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    ImageIndexType ind = it.GetIndex();

    bool isborder = false;
    for (unsigned int d = 0; d < PointDim; d++)
      if (ind[d] <= m_GridPadding || ind[d] >= ((long) m_GridSize[d] - m_GridPadding))
        isborder = true;
    if (isborder)
      continue;

    ImagePointType p;
    kernelImg->TransformIndexToPhysicalPoint(ind, p);

    VectorType x(PointDim, 0.0);
    for (unsigned int d = 0; d < PointDim; d++)
      x[d] = p[d];

    /*
  VectorType dvec = x - center;
  if (dvec.squared_magnitude() < nearThres)
    continue;
     */

    ScalarType k = this->EvaluateKernel(x, center);

    it.Set(k);
  }


  // Shift and flip halfplanes
  typedef itk::FFTShiftImageFilter<ImageType, ImageType>
      ShifterType;
  typename ShifterType::Pointer shiftf = ShifterType::New();
  shiftf->SetInput(kernelImg);
  shiftf->SetInverse(true);
  shiftf->SetNumberOfThreads(1);
  shiftf->Update();
  kernelImg = shiftf->GetOutput();


#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#endif

  //P3MKernel::m_FFTMutex.Lock();

  kernelImg->ReleaseDataFlagOff();

  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
  forwardFFT->SetInput(kernelImg);
  forwardFFT->SetNumberOfThreads(1);
  forwardFFT->Update();


  //P3MKernel::m_FFTMutex.Unlock();

  ComplexImagePointer kernelFFTImage = forwardFFT->GetOutput();

  P3MKernel::m_CacheMutex.Lock();
  m_CacheFFTKernelWidths.push_back(kernelWidth);
  m_CacheFFTKernelSizes.push_back(m_GridSize);
  m_CacheFFTKernelImages.push_back(kernelFFTImage);
  P3MKernel::m_CacheMutex.Unlock();

  return kernelFFTImage;
}

template<class ScalarType, unsigned int PointDim>
std::vector<typename P3MKernel<ScalarType, PointDim>::ComplexImagePointer>
P3MKernel<ScalarType, PointDim>
::BuildFFTGradientKernels() {
  ScalarType kernelWidth = this->GetKernelWidth();

  P3MKernel::m_CacheMutex.Lock();

  unsigned int numCache = m_CacheFFTGradientKernelWidths.size();

  unsigned int cacheIndex = numCache + 1;
  for (unsigned int i = 0; i < numCache; i++)
    if (m_CacheFFTGradientKernelWidths[i] == kernelWidth) {
      bool sameSize = true;
      ImageSizeType size_i = m_CacheFFTGradientKernelSizes[i];
      for (unsigned int d = 0; d < PointDim; d++)
        if (size_i[d] != m_GridSize[d]) {
          sameSize = false;
          break;
        }
      if (sameSize) {
        cacheIndex = i;
        break;
      }
    }

  P3MKernel::m_CacheMutex.Unlock();

  if (cacheIndex < numCache) {
    return m_CacheFFTGradientKernelImages[cacheIndex];
  }

  //std::cout << "Building grad kernel with width = " << kernelWidth << " and size = " << m_GridSize << std::endl;

  std::vector<ImagePointer> kernList;

  for (unsigned int dim = 0; dim < PointDim; dim++) {
    ImagePointer kernelImg = ImageType::New();
    kernelImg->SetRegions(m_MeshList[0]->GetLargestPossibleRegion());
    kernelImg->Allocate();
    kernelImg->CopyInformation(m_MeshList[0]);
    kernelImg->FillBuffer(0);

    kernList.push_back(kernelImg);
  }

  VectorType center(PointDim, 0.0);

//	ScalarType sumK = 1e-20;

  for (unsigned int d = 0; d < PointDim; d++) {
    center[d] =
        m_GridOrigin[d] + 0.5f * m_GridSize[d] * m_GridSpacing[d];
  }

  ImageIteratorType it(
      kernList[0], kernList[0]->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    ImageIndexType ind = it.GetIndex();

    bool isborder = false;
    for (unsigned int d = 0; d < PointDim; d++)
      if (ind[d] <= m_GridPadding || ind[d] >= ((long) m_GridSize[d] - m_GridPadding))
        isborder = true;
    if (isborder)
      continue;

    ImagePointType p;
    kernList[0]->TransformIndexToPhysicalPoint(ind, p);

    VectorType x(PointDim, 0.0);
    for (unsigned int d = 0; d < PointDim; d++)
      x[d] = p[d];

    VectorType g = this->EvaluateKernelGradient(x, center);

    for (unsigned int dim = 0; dim < PointDim; dim++)
      kernList[dim]->SetPixel(ind, g[dim]);

    // ScalarType k = this->EvaluateKernel(x, center);
    //
    // sumK += k;
  }

  // Shift and flip halfplanes
  for (unsigned int dim = 0; dim < PointDim; dim++) {
    typedef itk::FFTShiftImageFilter<ImageType, ImageType>
        ShifterType;
    typename ShifterType::Pointer shiftf = ShifterType::New();
    shiftf->SetInput(kernList[dim]);
    shiftf->SetInverse(true);
    shiftf->Update();
    kernList[dim] = shiftf->GetOutput();
  }

#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#endif

  std::vector<ComplexImagePointer> gradFFTKernels;

  for (unsigned int dim = 0; dim < PointDim; dim++) {
    //P3MKernel::m_FFTMutex.Lock();

    kernList[dim]->ReleaseDataFlagOff();

    typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
    forwardFFT->SetInput(kernList[dim]);
    forwardFFT->SetNumberOfThreads(1);
    forwardFFT->Update();

    //P3MKernel::m_FFTMutex.Unlock();

    gradFFTKernels.push_back(forwardFFT->GetOutput());
  } // for dim

  P3MKernel::m_CacheMutex.Lock();
  m_CacheFFTGradientKernelWidths.push_back(kernelWidth);
  m_CacheFFTGradientKernelSizes.push_back(m_GridSize);
  m_CacheFFTGradientKernelImages.push_back(gradFFTKernels);
  P3MKernel::m_CacheMutex.Unlock();

  return gradFFTKernels;
}

template<class ScalarType, unsigned int PointDim>
std::vector<typename P3MKernel<ScalarType, PointDim>::ComplexImagePointer>
P3MKernel<ScalarType, PointDim>
::BuildFFTHessianKernels() {
  ScalarType kernelWidth = this->GetKernelWidth();

  P3MKernel::m_CacheMutex.Lock();

  unsigned int numCache = m_CacheFFTHessianKernelWidths.size();

  unsigned int cacheIndex = numCache + 1;
  for (unsigned int i = 0; i < numCache; i++)
    if (m_CacheFFTHessianKernelWidths[i] == kernelWidth) {
      bool sameSize = true;
      ImageSizeType size_i = m_CacheFFTHessianKernelSizes[i];
      for (unsigned int d = 0; d < PointDim; d++)
        if (size_i[d] != m_GridSize[d]) {
          sameSize = false;
          break;
        }
      if (sameSize) {
        cacheIndex = i;
        break;
      }
    }

  P3MKernel::m_CacheMutex.Unlock();

  if (cacheIndex < numCache) {
    return m_CacheFFTHessianKernelImages[cacheIndex];
  }

  //std::cout << "Building hess kernel with width = " << kernelWidth << " and size = " << m_GridSize << std::endl;

  std::vector<ImagePointer> kernList;

  for (unsigned int r = 0; r < PointDim; r++)
    for (unsigned int c = r; c < PointDim; c++) {
      ImagePointer kernelImg = ImageType::New();
      kernelImg->SetRegions(m_MeshList[0]->GetLargestPossibleRegion());
      kernelImg->Allocate();
      kernelImg->CopyInformation(m_MeshList[0]);
      kernelImg->FillBuffer(0);

      kernList.push_back(kernelImg);
    }

  VectorType center(PointDim, 0.0);

//	ScalarType sumK = 1e-20;

  for (unsigned int d = 0; d < PointDim; d++) {
    center[d] =
        m_GridOrigin[d] + 0.5f * m_GridSize[d] * m_GridSpacing[d];
  }

  ImageIteratorType it(
      kernList[0], kernList[0]->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    ImageIndexType ind = it.GetIndex();

    bool isborder = false;
    for (unsigned int d = 0; d < PointDim; d++)
      if (ind[d] <= m_GridPadding || ind[d] >= ((long) m_GridSize[d] - m_GridPadding))
        isborder = true;
    if (isborder)
      continue;

    ImagePointType p;
    kernList[0]->TransformIndexToPhysicalPoint(ind, p);

    VectorType x(PointDim, 0.0);
    for (unsigned int d = 0; d < PointDim; d++)
      x[d] = p[d];

    MatrixType H = this->EvaluateKernelHessian(x, center);

    unsigned int ikernel = 0;
    for (unsigned int r = 0; r < PointDim; r++)
      for (unsigned int c = r; c < PointDim; c++)
        kernList[ikernel++]->SetPixel(ind, H(r, c));

    // ScalarType k = this->EvaluateKernel(x, center);
    //
    // sumK += k;
  }

  // Shift and flip halfplanes
  for (unsigned int j = 0; j < kernList.size(); j++) {
    typedef itk::FFTShiftImageFilter<ImageType, ImageType>
        ShifterType;
    typename ShifterType::Pointer shiftf = ShifterType::New();
    shiftf->SetInput(kernList[j]);
    shiftf->SetInverse(true);
    shiftf->Update();
    kernList[j] = shiftf->GetOutput();
  }
#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
#endif

#endif

  std::vector<ComplexImagePointer> hessFFTKernels;

  for (unsigned int j = 0; j < kernList.size(); j++) {
    //P3MKernel::m_FFTMutex.Lock();

    kernList[j]->ReleaseDataFlagOff();

    typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
    forwardFFT->SetInput(kernList[j]);
    forwardFFT->SetNumberOfThreads(1);
    forwardFFT->Update();

    //P3MKernel::m_FFTMutex.Unlock();

    hessFFTKernels.push_back(forwardFFT->GetOutput());
  } // for dim

  P3MKernel::m_CacheMutex.Lock();
  m_CacheFFTHessianKernelWidths.push_back(kernelWidth);
  m_CacheFFTHessianKernelSizes.push_back(m_GridSize);
  m_CacheFFTHessianKernelImages.push_back(hessFFTKernels);
  P3MKernel::m_CacheMutex.Unlock();

  return hessFFTKernels;
}

template<class ScalarType, unsigned int PointDim>
typename P3MKernel<ScalarType, PointDim>::ImagePointer
P3MKernel<ScalarType, PointDim>
::ApplyKernelFFT(ComplexImageType *kernelImg, ImageType *img) {
  if (kernelImg == 0) {

    std::cout << "no kernelImg!... implementation doubtful..." << std::endl;

    ScalarType h = this->GetKernelWidth();

    //typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>
    typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>
        GaussFilterType;
    typename GaussFilterType::Pointer gaussf = GaussFilterType::New();
    gaussf->SetInput(img);
    //gaussf->SetVariance(h*h);
    gaussf->SetSigma(h);
    gaussf->SetNormalizeAcrossScale(false);
    gaussf->Update();

    ImagePointer outimg = gaussf->GetOutput();

    ScalarType Z = pow(2.0f * vnl_math::pi * h, (ScalarType) PointDim / 2);

    ImageIteratorType it(outimg, outimg->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      it.Set(it.Get() * Z);

    return outimg;
  }

#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::FFTWComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::FFTWInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::VnlInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#endif

  // TODO:
  // fftw execute is thread safe, but planner is not
  // for thread safety and speed may need to write own ITK fftw wrapper

  //P3MKernel::m_FFTMutex.Lock();

  img->ReleaseDataFlagOff();

  typename ForwardFFTType::Pointer forwardFFT = ForwardFFTType::New();
  forwardFFT->SetInput(img);
  forwardFFT->SetNumberOfThreads(1);
  forwardFFT->Update();

  //P3MKernel::m_FFTMutex.Unlock();

  ComplexImagePointer fftImg = forwardFFT->GetOutput();
  fftImg->ReleaseDataFlagOff();

  typedef itk::ImageRegionIteratorWithIndex<ComplexImageType>
      ComplexImageIteratorType;

  ComplexImageIteratorType f_it(fftImg, fftImg->GetLargestPossibleRegion());
  for (f_it.GoToBegin(); !f_it.IsAtEnd(); ++f_it) {
    ImageIndexType ind = f_it.GetIndex();
    f_it.Set(f_it.Get() * kernelImg->GetPixel(ind));
  }

  //P3MKernel::m_FFTMutex.Lock();

  typename InverseFFTType::Pointer invFFT = InverseFFTType::New();
  invFFT->SetInput(fftImg);
  invFFT->SetNumberOfThreads(1);
  invFFT->Update();

  //P3MKernel::m_FFTMutex.Unlock();

  ImagePointer outimg = invFFT->GetOutput();

  return outimg;
}

template<class ScalarType, unsigned int PointDim>
typename P3MKernel<ScalarType, PointDim>::ImagePointer
P3MKernel<ScalarType, PointDim>
::ApplyKernelFFT(ComplexImageType *kernelImg, ComplexImageType *img) {
  if (kernelImg == 0)
    throw std::runtime_error("should give kernelImg");

#if defined(USE_FFTWF)

#if ITK_VERSION_MAJOR < 4
  typedef itk::FFTWRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::FFTWComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::FFTWForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::FFTWInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#else

#if ITK_VERSION_MAJOR < 4
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<float, PointDim> ForwardFFTType;
  typedef itk::VnlFFTComplexConjugateToRealImageFilter<float, PointDim> InverseFFTType;
#else
  typedef itk::VnlForwardFFTImageFilter<ImageType> ForwardFFTType;
  typedef itk::VnlInverseFFTImageFilter<ComplexImageType> InverseFFTType;
#endif

#endif

  typedef itk::ImageDuplicator<ComplexImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(img);
  duplicator->Update();
  typename ComplexImageType::Pointer AuxImg = duplicator->GetOutput();
  AuxImg->ReleaseDataFlagOff();

  typedef itk::ImageRegionIteratorWithIndex<ComplexImageType>
      ComplexImageIteratorType;

  ComplexImageIteratorType f_it(AuxImg, AuxImg->GetLargestPossibleRegion());
  for (f_it.GoToBegin(); !f_it.IsAtEnd(); ++f_it) {
    ImageIndexType ind = f_it.GetIndex();
    f_it.Set(f_it.Get() * kernelImg->GetPixel(ind));
  }

  // TODO:
  // fftw execute is thread safe, but planner is not
  // for thread safety and speed may need to write own ITK fftw wrapper

  //P3MKernel::m_FFTMutex.Lock();

  typename InverseFFTType::Pointer invFFT = InverseFFTType::New();
  invFFT->SetInput(AuxImg);
  invFFT->SetNumberOfThreads(1);
  invFFT->Update();

  //P3MKernel::m_FFTMutex.Unlock();

  ImagePointer outimg = invFFT->GetOutput();

  return outimg;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
P3MKernel<ScalarType, PointDim>
::Convolve(const MatrixType &X) {
//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//    auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//    m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  MatrixType V(X.rows(), weightDim, 0.0);

  // ScalarType nearThres = 0;
  // for (unsigned int d = 0; d < PointDim; d++)
  //   nearThres += m_GridSpacing[d]*m_GridSpacing[d];
  // nearThres *= m_NearThresholdScale*m_NearThresholdScale;

//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  // Convolve weights splatted on mesh with kernel
  std::vector<ImagePointer> img(weightDim);
  for (unsigned int k = 0; k < weightDim; k++)
    img[k] = this->ApplyKernelFFT(m_FFTKernel, m_MeshListFFT[k]);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//    auto da = std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count();
//    m_ApplyKernelFftTime += da;

  for (unsigned int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);

    VectorType vi = this->Interpolate(xi, img);

#if DO_NEAR_FIELD
    //long id = this->_getPointID(xi);

    std::vector<ScalarType> weights;
    std::vector<ImageIndexType> gridIndices;
    this->_getInterpolationWeightsAndGridPoints(
      weights, gridIndices, xi);

    std::vector<long> idlist;
    for (unsigned int j = 0; j < gridIndices.size(); j++)
      idlist.push_back(this->_getPointID(gridIndices[j]));

    for (unsigned int j = 0; j < m_SourceIDs.size(); j++)
    {
#if DO_SORT_ID
      if (m_SourceIDs[j] > id)
        break;
#endif

      //if (m_SourceIDs[j] != id)
      //  continue;
      bool isneigh = false;
      for (unsigned int k = 0; k < idlist.size(); k++)
        if (idlist[k] == m_SourceIDs[j])
        {
          isneigh = true;
          break;
        }
      if (!isneigh)
        continue;

      VectorType yj = Superclass::m_Sources.get_row(j);

      VectorType dvec = xi - yj;

      if (dvec.squared_magnitude() < nearThres)
      {
        ScalarType sij = this->EvaluateKernel(xi, yj);
        for (unsigned int k = 0; k < Superclass::m_Weights.columns(); k++)
          vi[k] += sij * Superclass::m_Weights(j, k);
      }
    }

#endif

    V.set_row(i, vi);
  }

//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//	auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
//	m_ConvolveTime += d;
//
//    auto di = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t4).count();
//    m_InterpolateTime += di;

  return V;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
P3MKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, const MatrixType &alpha) {
  std::vector<MatrixType> gradMom = this->ConvolveGradient(X);

  MatrixType result(X.rows(), PointDim, 0);
  for (unsigned int j = 0; j < X.rows(); j++)
    result.set_row(j, gradMom[j].transpose() * alpha.get_row(j));

  return result;
}

template<class ScalarType, unsigned int PointDim>
std::vector<MatrixType>
P3MKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X) {
//     /// For profiling.
//     std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

//     /// For profiling.
//     std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//     auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//     m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  // ScalarType nearThres = 0;
  // for (unsigned int d = 0; d < PointDim; d++)
  //   nearThres += m_GridSpacing[d]*m_GridSpacing[d];
  // nearThres *= m_NearThresholdScale*m_NearThresholdScale;

  std::vector<MatrixType> gradK(
      X.rows(), MatrixType(weightDim, PointDim, 0.0));

//     /// For profiling.
//     std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  std::vector<ImagePointer> img(weightDim * PointDim);
  for (unsigned int k = 0; k < weightDim; k++)
    for (unsigned int p = 0; p < PointDim; p++)
      img[p + PointDim * k] = this->ApplyKernelFFT(m_FFTGradientKernels[p], m_MeshListFFT[k]);

//     /// For profiling.
//     std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//     auto da = std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count();
//     m_ApplyKernelFftTime += da;

  for (unsigned int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);

    VectorType wi = this->Interpolate(xi, img);

    MatrixType &G = gradK[i];
    for (unsigned int k = 0; k < weightDim; k++)
      for (unsigned int p = 0; p < PointDim; p++)
        G(k, p) = wi[p + PointDim * k];

  }

//     /// For profiling.
//     std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//     auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
//     m_ConvolveGradientTime += d;
//
//     auto di = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t4).count();
//     m_InterpolateTime += di;

  return gradK;

}

template<class ScalarType, unsigned int PointDim>
MatrixType
P3MKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, unsigned int dim) {
//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//    auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//    m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");
  if ((dim < 0) || (dim >= PointDim))
    throw std::runtime_error("dimension parameter out of bounds");

  MatrixType gradK(X.rows(), weightDim, 0.0);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  std::vector<ImagePointer> img(weightDim);
  for (unsigned int k = 0; k < weightDim; k++)
    img[k] = this->ApplyKernelFFT(m_FFTGradientKernels[dim], m_MeshListFFT[k]);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//    auto da = std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count();
//    m_ApplyKernelFftTime += da;

  for (unsigned int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);

    VectorType wi = this->Interpolate(xi, img);

    gradK.set_row(i, wi);
  }

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//    auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
//    m_ConvolveGradientTime += d;
//
//    auto di = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t4).count();
//    m_InterpolateTime += di;

  return gradK;
}

template<class ScalarType, unsigned int PointDim>
VectorType
P3MKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, unsigned int k, unsigned int dp) {
//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//    auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//    m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  if (k >= weightDim)
    throw std::runtime_error("Invalid weight index");
  if (dp >= PointDim)
    throw std::runtime_error("Invalid derivative direction");

  VectorType gradK(X.rows(), 0.0);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//    auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
//    m_ConvolveGradientTime += d;

  return gradK;
}

template<class ScalarType, unsigned int PointDim>
std::vector<std::vector<MatrixType> >
P3MKernel<ScalarType, PointDim>
::ConvolveHessian(const MatrixType &X) {
//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

  if (!m_HessianUpdated)
    this->UpdateHessianGrids();

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//    auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//    m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  std::vector<std::vector<MatrixType> > hessK;
  for (unsigned int i = 0; i < X.rows(); i++) {
    MatrixType M(PointDim, PointDim, 0.0);
    std::vector<MatrixType> H(weightDim, M);
    hessK.push_back(H);
  }

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  int ptsD2 = PointDim * (PointDim + 1) / 2;
  std::vector<ImagePointer> img(weightDim * ptsD2);
  for (unsigned int p = 0; p < ptsD2; p++)
    for (unsigned int k = 0; k < weightDim; k++)
      img[k + weightDim * p] = this->ApplyKernelFFT(m_FFTHessianKernels[p], m_MeshListFFT[k]);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//    auto da = std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count();
//    m_ApplyKernelFftTime += da;

  for (unsigned int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);
    VectorType w = this->Interpolate(xi, img);

    for (unsigned int k = 0; k < weightDim; k++) {
      MatrixType Hik(PointDim, PointDim, 0.0);
      unsigned int idx = 0;
      for (unsigned int r = 0; r < PointDim; r++)
        for (unsigned int c = r; c < PointDim; c++) {
          Hik(r, c) = w[k + weightDim * idx];
          Hik(c, r) = Hik(r, c);
          idx++;
        }

      hessK[i][k] = Hik;
    }
  }

  return hessK;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
P3MKernel<ScalarType, PointDim>
::ConvolveHessian(const MatrixType &X, unsigned int row, unsigned int col) {
//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

  if (!m_HessianUpdated)
    this->UpdateHessianGrids();

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//    auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//    m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");
  if (row >= PointDim || col >= PointDim)
    throw std::runtime_error("dimension index out of bounds");

  MatrixType hessK(X.rows(), weightDim, 0.0);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  int index =
      (row <= col) ? (col + PointDim * row - row * (row + 1) / 2) : (row + PointDim * col - col * (col + 1) / 2);
  std::vector<ImagePointer> img(weightDim);
  for (unsigned int k = 0; k < weightDim; k++)
    img[k] = this->ApplyKernelFFT(m_FFTHessianKernels[index], m_MeshListFFT[k]);

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
//    auto da = std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count();
//    m_ApplyKernelFftTime += da;

  for (int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);
    VectorType wi = this->Interpolate(xi, img);

    hessK.set_row(i, wi);
  }

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//    auto d = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
//    m_ConvolveHessianTime += d;
//
//    auto di = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t4).count();
//    m_InterpolateTime += di;

  return hessK;
}

template<class ScalarType, unsigned int PointDim>
VectorType
P3MKernel<ScalarType, PointDim>
::ConvolveHessian(const MatrixType &X, unsigned int k,
                  unsigned int dp, unsigned int dq) {
//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  if (this->IsModified()) {
    this->UpdateGrids();
    this->UnsetModified();
  }

  if (!m_HessianUpdated)
    this->UpdateHessianGrids();

//    /// For profiling.
//    std::chrono::high_resolution_clock::time_point tg = std::chrono::high_resolution_clock::now();
//    auto dg = std::chrono::duration_cast<std::chrono::milliseconds>(tg-t1).count();
//    m_UpdateGridTime += dg;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  unsigned int sourceDim = Y.columns();
  unsigned int weightDim = W.columns();

  if (sourceDim != PointDim)
    throw std::runtime_error("Can only handle certain dimension");
  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  if (k >= weightDim)
    throw std::runtime_error("Invalid weight index");
  if (dp >= PointDim || dq >= PointDim)
    throw std::runtime_error("Invalid derivative direction");

  VectorType hessK(X.rows(), 0.0);

  return hessK;
}

template class P3MKernel<double, 2>;
template class P3MKernel<double, 3>;
