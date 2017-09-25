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

#include "CUDAExactKernel.h"
#include "SimpleTimer.h"
#include "../../lib/cuda_convolutions/GpuConv1D.h"
//#include "../../lib/cuda_convolutions/GpuConv2D.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int PointDim>
CUDAExactKernel<ScalarType, PointDim>
::CUDAExactKernel(const CUDAExactKernel &o) {
  Superclass::m_Sources = o.m_Sources;
  Superclass::m_Weights = o.m_Weights;
  this->SetKernelWidth(o.GetKernelWidth());

  if (o.IsModified())
    this->SetModified();
  else
    this->UnsetModified();
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int PointDim>
MatrixType
CUDAExactKernel<ScalarType, PointDim>
::Convolve(const MatrixType &X) {

  VectorType X_vec = X.vectorise_row_wise();
  VectorType S_vec = this->GetSources().vectorise_row_wise();
  VectorType W_vec = this->GetWeights().vectorise_row_wise();

  if (this->GetSources().rows() != this->GetWeights().rows())
    throw std::runtime_error("Sources and weights count mismatch");

  int DimVect = this->GetWeights().columns();
  ScalarType *gammap = new ScalarType[X.rows() * DimVect];

  if (DimVect == PointDim) {
    GaussGpuEvalConv1D<ScalarType, PointDim, PointDim>(this->GetKernelWidth(),
                                                       X_vec.memptr(),
                                                       S_vec.memptr(),
                                                       W_vec.memptr(),
                                                       gammap,
                                                       X.rows(),
                                                       this->GetSources().rows());
  } else if (DimVect == 2 * PointDim) {
    GaussGpuEvalConv1D<ScalarType, PointDim, (2 * PointDim)>(this->GetKernelWidth(),
                                                             X_vec.memptr(),
                                                             S_vec.memptr(),
                                                             W_vec.memptr(),
                                                             gammap,
                                                             X.rows(),
                                                             this->GetSources().rows());
  } else if (DimVect == PointDim * (PointDim + 1) / 2) {
    GaussGpuEvalConv1D<ScalarType, PointDim, PointDim * (PointDim + 1) / 2>(this->GetKernelWidth(),
                                                                            X_vec.memptr(),
                                                                            S_vec.memptr(),
                                                                            W_vec.memptr(),
                                                                            gammap,
                                                                            X.rows(),
                                                                            this->GetSources().rows());
  } else {
    throw std::runtime_error("In CUDAExactKernel::Convolve(X) - Invalid number of columns of beta !");
  }


  MatrixType gamma(gammap, DimVect, X.rows());
  //	return gamma.transpose();
  MatrixType transpose = gamma.transpose();
  delete[] gammap;
  return transpose;
}


//template<class ScalarType, unsigned int PointDim>
//MatrixType
//ExactKernel<ScalarType, PointDim>
//::ConvolveImageFast(const MatrixType &X, const ImageTypePointer image) {
//  typedef itk::Point<ScalarType, PointDim> PointType;
//  typedef itk::Index<PointDim> IndexType;
//  typedef itk::ImageRegion<PointDim> RegionType;
//  typedef itk::Size<PointDim> SizeType;
//
//  MatrixType &Y = this->GetSources(); // Control points.
//  MatrixType &W = this->GetWeights(); // Momenta.
//
//
//  if (Y.rows() != W.rows()) { throw std::runtime_error("Sources and weights count mismatch"); }
//
//  const SizeType sizeImage = image->GetLargestPossibleRegion().GetSize();
//  SizeType sizeRegion;
//  sizeRegion.Fill(6 * Superclass::m_KernelWidth);
//
//  MatrixType V(X.rows(), PointDim, 0.0);
//
//  unsigned int nbBatch = 20;
//  unsigned int currentCpInBatch = 0;
//
//  std::vector<unsigned int> indicesControlPoints;
//  std::vectr
//
//  for (unsigned int control_point_index = 0; control_point_index < Y.rows(); ++control_point_index) {
//    auto controlPoint = Y.get_row(control_point_index);
//
//    PointType regionReferencePoint;
//    for (unsigned int d = 0; d < PointDim; ++d) {
//      regionReferencePoint[d] = controlPoint[d] - 3 * Superclass::m_KernelWidth;
//    }
//
//    IndexType regionReferenceIndex;
//    image->TransformPhysicalPointToIndex(regionReferencePoint, regionReferenceIndex);
//
//    RegionType region(regionReferenceIndex, sizeRegion);
//    region.Crop(image->GetLargestPossibleRegion());
//
//    itk::ImageRegionConstIteratorWithIndex<ImageType> regionIterator(image, region);
//
//    while (!regionIterator.IsAtEnd()) {
//      IndexType currentIndex = regionIterator.GetIndex();
//      unsigned int pixel_index = currentIndex[0] + sizeImage[0] * currentIndex[1];
//      if (PointDim == 3) { pixel_index += sizeImage[0] * sizeImage[1] * currentIndex[2]; }
//
//      ScalarType Kij = this->EvaluateKernel(X, Y, pixel_index, control_point_index);
//      for (unsigned int d = 0; d < PointDim; ++d) { V(pixel_index, d) += W(control_point_index, d) * Kij; }
//
//      ++regionIterator;
//    }
//  }
//  return V;
//};
















template<class ScalarType, unsigned int PointDim>
MatrixType
CUDAExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, const MatrixType &alpha) {
  int DimVect = this->GetWeights().columns();
  ScalarType *gammap = new ScalarType[X.rows() * DimVect];

  VectorType alpha_vec = alpha.vectorise_row_wise();
  VectorType X_vec = X.vectorise_row_wise();
  VectorType S_vec = this->GetSources().vectorise_row_wise();
  VectorType W_vec = this->GetWeights().vectorise_row_wise();

  if (DimVect == PointDim) {
    GaussGpuGrad1Conv1D<ScalarType, PointDim, PointDim>(this->GetKernelWidth(),
                                                        alpha_vec.memptr(),
                                                        X_vec.memptr(),
                                                        S_vec.memptr(),
                                                        W_vec.memptr(),
                                                        gammap,
                                                        X.rows(),
                                                        this->GetSources().rows());
  } else {
    throw std::runtime_error(
        "In CUDAExactKernel::ConvolveGradient(X, alpha) - Problem with the number of columns of beta !");
  }
  MatrixType gamma(gammap, DimVect, X.rows());
  MatrixType transpose = gamma.transpose();
  delete[] gammap;
  return transpose;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
CUDAExactKernel<ScalarType, PointDim>
::ConvolveSpecialHessian(const MatrixType &xi) {
  if (this->GetSources().rows() != this->GetWeights().rows())
    throw std::runtime_error("Sources and weights count mismatch");
  if (this->GetSources().rows() != xi.rows())
    throw std::runtime_error("Y and Xi count mismatch");

  int DimVect = this->GetWeights().columns();
  ScalarType *gamma_p = new ScalarType[xi.rows() * DimVect];

  VectorType S_vec = this->GetSources().vectorise_row_wise();
  VectorType W_vec = this->GetWeights().vectorise_row_wise();
  VectorType xi_vec = xi.vectorise_row_wise();

  if (DimVect == PointDim) {
    GaussGpuGradDiffConv1D<ScalarType, PointDim, PointDim>(this->GetKernelWidth(),
                                                           S_vec.memptr(),
                                                           W_vec.memptr(),
                                                           xi_vec.memptr(),
                                                           gamma_p,
                                                           this->GetSources().rows());
  } else {
    throw std::runtime_error("In CUDAExactKernel::ConvolveSpecialHessian(xi) - Invalid number of columns of beta !");
  }
  MatrixType gamma(gamma_p, DimVect, xi.rows());
  auto ret = (-0.5f * gamma.transpose());
  delete[] gamma_p;
  return ret;
}

template<class ScalarType, unsigned int PointDim>
MatrixType
CUDAExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X, unsigned int dim) {

  VectorType X_vec = X.vectorise_row_wise();
  VectorType W_vec = this->GetWeights().vectorise_row_wise();

  int weightDim = this->GetWeights().columns();

  ScalarType *gamma_p = new ScalarType[X.rows() * weightDim];

  if (X.rows() != this->GetWeights().rows())
    throw std::runtime_error("Sources and weights count mismatch");
  if (dim >= PointDim || dim < 0)
    throw std::runtime_error("dimension index out of bounds");

  GaussGpuGradConv1D<ScalarType, PointDim, PointDim>(this->GetKernelWidth(),
                                                     X_vec.memptr(),
                                                     W_vec.memptr(),
                                                     dim,
                                                     gamma_p,
                                                     X.rows()
  );

  MatrixType gamma(gamma_p, weightDim, X.rows());
  MatrixType transpose = gamma.transpose();

  delete[] gamma_p;
  return transpose;

}

template<class ScalarType, unsigned int PointDim>
std::vector<MatrixType>
CUDAExactKernel<ScalarType, PointDim>
::ConvolveGradient(const MatrixType &X) {

  VectorType X_vec = X.vectorise_row_wise();
  VectorType Y_vec = this->GetSources().vectorise_row_wise();
  VectorType W_vec = this->GetWeights().vectorise_row_wise();

  if (this->GetSources().rows() != this->GetWeights().rows())
    throw std::runtime_error("Sources and weights count mismatch");

  unsigned int DimVect = this->GetWeights().columns();

  ScalarType *gamma_p = new ScalarType[X.rows() * DimVect * PointDim];

  if (DimVect == PointDim) {
    GaussGpuGradConv_varlin_1D<ScalarType, PointDim, PointDim>(this->GetKernelWidth(),
                                                               X_vec.memptr(),
                                                               Y_vec.memptr(),
                                                               W_vec.memptr(),
                                                               gamma_p,
                                                               X.rows(),
                                                               this->GetSources().rows());
  } else if (DimVect == 2 * PointDim) {
    GaussGpuGradConv_varlin_1D<ScalarType, PointDim, (2 * PointDim)>(this->GetKernelWidth(),
                                                                     X_vec.memptr(),
                                                                     Y_vec.memptr(),
                                                                     W_vec.memptr(),
                                                                     gamma_p,
                                                                     X.rows(),
                                                                     this->GetSources().rows());
  } else if (DimVect == PointDim * (PointDim + 1) / 2) {
    GaussGpuGradConv_varlin_1D<ScalarType, PointDim, PointDim * (PointDim + 1) / 2>(this->GetKernelWidth(),
                                                                                    X_vec.memptr(),
                                                                                    Y_vec.memptr(),
                                                                                    W_vec.memptr(),
                                                                                    gamma_p,
                                                                                    X.rows(),
                                                                                    this->GetSources().rows());
  } else {
    throw std::runtime_error("In CUDAExactKernel::ConvolveGradient(X) - Invalid number of columns of beta !");
  }

  std::vector<MatrixType> gradK;

  for (unsigned int i = 0; i < X.rows(); i++) {
    MatrixType Gi(DimVect, PointDim, 0.0);

    for (unsigned int k = 0; k < DimVect; k++) {
      for (unsigned int l = 0; l < PointDim; l++) {
        Gi(k, l) = gamma_p[i * DimVect * PointDim + k * PointDim + l];
      }
    }
    gradK.push_back(Gi);
  }

  delete[] gamma_p;
  return gradK;
}

template class CUDAExactKernel<double, 2>;
template class CUDAExactKernel<double, 3>;

