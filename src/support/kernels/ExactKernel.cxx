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

#include "ExactKernel.h"

#include <exception>
#include <stdexcept>
#include <itkImageRegionIterator.h>

#include "SimpleTimer.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int PointDim>
ExactKernel<ScalarType, PointDim>
::ExactKernel(const ExactKernel &o) {
  Superclass::m_Sources = o.m_Sources;
  Superclass::m_Weights = o.m_Weights;
  this->SetKernelWidth(o.GetKernelWidth());

  if (o.IsModified())
    this->SetModified();
  else
    this->UnsetModified();

}

template<class ScalarType, unsigned int PointDim>
ExactKernel<ScalarType, PointDim> *
ExactKernel<ScalarType, PointDim>
::Clone() const {
  return new ExactKernel(*this);
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::Convolve(const MatrixType &X) {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  MatrixType V(X.rows(), weightDim, 0.0);

  for (unsigned int i = 0; i < X.rows(); i++) {
    for (unsigned int j = 0; j < Y.rows(); j++) {
      ScalarType Kij = this->EvaluateKernel(X, Y, i, j);
      for (unsigned int k = 0; k < weightDim; k++)
        V(i, k) += W(j, k) * Kij;
    }
  }

  return V;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ConvolveImageFast(const MatrixType &X, const ImageTypePointer image) {
  typedef itk::Point<ScalarType, PointDim> PointType;
  typedef itk::Index<PointDim> IndexType;
  typedef itk::ImageRegion<PointDim> RegionType;
  typedef itk::Size<PointDim> SizeType;

  MatrixType &Y = this->GetSources(); // Control points.
  MatrixType &W = this->GetWeights(); // Momenta.

  std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

  if (Y.rows() != W.rows()) { throw std::runtime_error("Sources and weights count mismatch"); }

  const SizeType sizeImage = image->GetLargestPossibleRegion().GetSize();

  SizeType sizeRegion;
  sizeRegion.Fill(6 * Superclass::m_KernelWidth);

  MatrixType V(X.rows(), PointDim, 0.0);

  for (unsigned int control_point_index = 0; control_point_index < Y.rows(); ++control_point_index) {
    auto controlPoint = Y.get_row(control_point_index);

    PointType regionReferencePoint;
    for (unsigned int d = 0; d < PointDim; ++d) {
      regionReferencePoint[d] = controlPoint[d] - 3 * Superclass::m_KernelWidth;
    }

    IndexType regionReferenceIndex;
    image->TransformPhysicalPointToIndex(regionReferencePoint, regionReferenceIndex);

    RegionType region(regionReferenceIndex, sizeRegion);
    region.Crop(image->GetLargestPossibleRegion());

    itk::ImageRegionConstIteratorWithIndex<ImageType> regionIterator(image, region);

    while (!regionIterator.IsAtEnd()) {
      IndexType currentIndex = regionIterator.GetIndex();
      unsigned int pixel_index = currentIndex[0] + sizeImage[0] * currentIndex[1];
      if (PointDim == 3) { pixel_index += sizeImage[0] * sizeImage[1] * currentIndex[2]; }

      ScalarType Kij = this->EvaluateKernel(X, Y, pixel_index, control_point_index);
      for (unsigned int d = 0; d < PointDim; ++d) { V(pixel_index, d) += W(control_point_index, d) * Kij; }

      ++regionIterator;
    }
  }

  return V;

};

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ComputeKernelMatrix(const MatrixType &Y) {
  const unsigned int N = Y.rows();
  MatrixType matKernel(N, N, 0.);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      matKernel(i, j) = this->EvaluateKernel(Y, Y, i, j);
    }
  }
  return matKernel;
};

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ComputeInverseMatrix(const MatrixType &M) {
  MatrixType invM;//TODO : check the conditionning of the matrix and look for pseudo-inverse ?
  invM = inverse_sympd(M);
  return invM;
};

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ConvolveInverse(const MatrixType &X) {
  MatrixType &sources = this->GetSources();
  MatrixType &weights = this->GetWeights();
  MatrixType kernelMatrix(X.rows(), X.rows(), 0);
  kernelMatrix = this->ComputeKernelMatrix(sources);
  MatrixType V = solve(kernelMatrix, X);
  return V;
};

template<class ScalarType, unsigned int PointDim>
std::vector<MatrixType>
ExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X) {

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  std::vector<MatrixType> gradK;

  for (unsigned int i = 0; i < X.rows(); i++) {
    // 	VectorType xi = X.get_row(i);

    MatrixType Gi(weightDim, PointDim, 0.0);

    for (unsigned int j = 0; j < Y.rows(); j++) {
      // VectorType yj = Y.get_row(j);

      VectorType g = this->EvaluateKernelGradient(X, Y, i, j);
      for (unsigned int k = 0; k < weightDim; k++) {
        ScalarType Wjk = W(j, k);
        for (unsigned int l = 0; l < PointDim; l++) {
          Gi(k, l) += g[l] * Wjk;
        }
      }
    }

    gradK.push_back(Gi);
  }

  return gradK;
}

template<class ScalarType, unsigned int PointDim>
std::vector<MatrixType>
ExactKernel<ScalarType, PointDim>
::ConvolveGradientImageFast(const MatrixType &X, const ImageTypePointer image) {

  std::cout << "In convolvegradientimagefast(x,image)" << std::endl;

  typedef itk::Point<ScalarType, PointDim> PointType;
  typedef itk::Index<PointDim> IndexType;
  typedef itk::ImageRegion<PointDim> RegionType;
  typedef itk::Size<PointDim> SizeType;

  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  std::vector<MatrixType> gradK(X.rows(), MatrixType(weightDim, PointDim, 0.0));

  const SizeType sizeImage = image->GetLargestPossibleRegion().GetSize();
  SizeType sizeRegion;
  sizeRegion.Fill(6 * Superclass::m_KernelWidth);

  for (unsigned int control_point_index = 0; control_point_index < Y.rows(); ++control_point_index) {
    auto controlPoint = Y.get_row(control_point_index);

    PointType regionReferencePoint;
    for (unsigned int d = 0; d < PointDim; ++d) {
      regionReferencePoint[d] = controlPoint[d] - 3 * Superclass::m_KernelWidth;
    }

    IndexType regionReferenceIndex;
    image->TransformPhysicalPointToIndex(regionReferencePoint, regionReferenceIndex);

    RegionType region(regionReferenceIndex, sizeRegion);
    region.Crop(image->GetLargestPossibleRegion());

    itk::ImageRegionConstIteratorWithIndex<ImageType> regionIterator(image, region);

    while (!regionIterator.IsAtEnd()) {
      IndexType currentIndex = regionIterator.GetIndex();
      unsigned int pixel_index = currentIndex[0] + sizeImage[0] * currentIndex[1];
      if (PointDim == 3) { pixel_index += sizeImage[0] * sizeImage[1] * currentIndex[2]; }

      VectorType g = this->EvaluateKernelGradient(X, Y, pixel_index, control_point_index);
      for (unsigned int k = 0; k < weightDim; ++k) {
        ScalarType Wjk = W(control_point_index, k);
        for (unsigned int l = 0; l < PointDim; ++l) { gradK[pixel_index](k, l) += g[l] * Wjk; }
      }
      ++regionIterator;
    }
  }
  return gradK;
}



template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, const MatrixType &alpha) {
  std::vector<MatrixType> convolveGradient = this->ConvolveGradient(X);

  MatrixType result(X.rows(), PointDim, 0);
  for (unsigned int j = 0; j < X.rows(); j++)
    result.set_row(j, convolveGradient[j].transpose() * alpha.get_row(j));

  return result;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ConvolveGradientImageFast(const MatrixType &X, const MatrixType &alpha, const ImageTypePointer img) {

  std::cout << "In convolvegradientimagefast(x,alpha,image)" << std::endl;


  std::vector<MatrixType> convolveGradient = this->ConvolveGradientImageFast(X, img);

  MatrixType result(X.rows(), PointDim, 0);
  for (unsigned int j = 0; j < X.rows(); j++)
    result.set_row(j, convolveGradient[j].transpose() * alpha.get_row(j));

  return result;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, unsigned int dim) {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");
  if (dim >= PointDim || dim < 0)
    throw std::runtime_error("dimension index out of bounds");

  unsigned int weightDim = W.columns();

  MatrixType gradK(X.rows(), weightDim, 0.0);

  for (unsigned int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);

    for (unsigned int j = 0; j < Y.rows(); j++) {
      VectorType yj = Y.get_row(j);
      VectorType wj = W.get_row(j);

      VectorType g = this->EvaluateKernelGradient(xi, yj);
      gradK.set_row(i, gradK.get_row(i) + g[dim] * wj);
    }
  }

  return gradK;
}

template<class ScalarType, unsigned int PointDim>
VectorType
ExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, unsigned int k, unsigned int dp) {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  if (k >= weightDim)
    throw std::runtime_error("Invalid weight index");

  if (dp >= PointDim)
    throw std::runtime_error("Invalid derivative direction");

  unsigned int numPoints = Y.rows();

  VectorType gradK(numPoints, 0);

  for (unsigned int i = 0; i < numPoints; i++) {
    VectorType xi = X.get_row(i);

    gradK[i] = 0;

    for (unsigned int j = 0; j < numPoints; j++) {
      VectorType yj = Y.get_row(j);
      VectorType wj = W.get_row(j);

      VectorType g = this->EvaluateKernelGradient(xi, yj);

      gradK[i] += wj[k] * g[dp];
    }
  }

  return gradK;
}

template<class ScalarType, unsigned int PointDim>
std::vector<std::vector<MatrixType> >
ExactKernel<ScalarType, PointDim>
::ConvolveHessian(const MatrixType &X) {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  std::vector<std::vector<MatrixType> > hessK;

  for (unsigned int i = 0; i < X.rows(); i++) {

    std::vector<MatrixType> Hi;
    for (unsigned int k = 0; k < weightDim; k++)
      Hi.push_back(MatrixType(PointDim, PointDim, 0.0));

    for (unsigned int j = 0; j < Y.rows(); j++) {
      VectorType wj = W.get_row(j);

      MatrixType H = this->EvaluateKernelHessian(X, Y, i, j);
      for (unsigned int k = 0; k < weightDim; k++)
        Hi[k] = Hi[k] + H * wj[k];
    }

    hessK.push_back(Hi);
  }

  return hessK;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
ExactKernel<ScalarType, PointDim>
::ConvolveHessian(const MatrixType &X, unsigned int row, unsigned int col) {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");
  if (row >= PointDim || col >= PointDim)
    throw std::runtime_error("Dimension index out of bounds");

  unsigned int weightDim = W.columns();

  MatrixType hessK(X.rows(), weightDim, 0.0);

  for (unsigned int i = 0; i < X.rows(); i++) {
    VectorType xi = X.get_row(i);

    for (unsigned int j = 0; j < Y.rows(); j++) {
      VectorType yj = Y.get_row(j);
      VectorType wj = W.get_row(j);

      MatrixType H = this->EvaluateKernelHessian(xi, yj);
      hessK.set_row(i, hessK.get_row(i) + wj * H(row, col));
    }
  }

  return hessK;
}

template<class ScalarType, unsigned int PointDim>
VectorType
ExactKernel<ScalarType, PointDim>
::ConvolveHessian(const MatrixType &X, unsigned int k,
                  unsigned int dp, unsigned int dq) {
  MatrixType &Y = this->GetSources();
  MatrixType &W = this->GetWeights();

  if (Y.rows() != W.rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int weightDim = W.columns();

  if (k >= weightDim)
    throw std::runtime_error("Invalid weight index");

  if (dp >= PointDim || dq >= PointDim)
    throw std::runtime_error("Invalid derivative direction");

  unsigned int numPoints = X.rows();

  VectorType hessK(numPoints, 0);

  for (unsigned int i = 0; i < numPoints; i++) {
    VectorType xi = X.get_row(i);

    hessK[i] = 0;

    for (unsigned int j = 0; j < Y.rows(); j++) {
      VectorType yj = Y.get_row(j);
      VectorType wj = W.get_row(j);

      MatrixType H = this->EvaluateKernelHessian(xi, yj);

      hessK[i] += wj[k] * H(dp, dq);

    }
  }

  return hessK;
}

template
class ExactKernel<double, 2>;
template
class ExactKernel<double, 3>;