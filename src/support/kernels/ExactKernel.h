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

#pragma once

#include <cmath>
#include <exception>
#include <stdexcept>

#include <itkImage.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include "AbstractKernel.h"
#include "MathFunctions.h"

/**
 *	\brief      An exact kernel.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    The ExactKernel class inherited from AbstractKernel implements the operations of
 *              convolution and evaluation of the kernel using an exact computation.
 */
template<class ScalarType, unsigned int PointDim>
class ExactKernel : public AbstractKernel<ScalarType, PointDim> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract kernel type.
  typedef AbstractKernel<ScalarType, PointDim> Superclass;
  /// ITK image type.
  typedef itk::Image<ScalarType, PointDim> ImageType;
  /// ITK image pointer type.
  typedef typename ImageType::Pointer ImageTypePointer;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ExactKernel() : Superclass() {}
  /// Copy constructor.
  ExactKernel(const ExactKernel &o);
  /// See AbstractKernel::AbstractKernel(const MatrixType& X, double h).
  ExactKernel(const MatrixType &X, double h) : Superclass(X, h) {}
  /// See AbstractKernel::AbstractKernel(const MatrixType& X, const MatrixType& W, double h).
  ExactKernel(const MatrixType &X, const MatrixType &W, double h) : Superclass(X, W, h) {}

  virtual ExactKernel *Clone() const;

  virtual ~ExactKernel() {}



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual ScalarType EvaluateKernel(const VectorType &x, const VectorType &y) {
    // VectorType dvec = x - y;
    ScalarType distsq = (x - y).squared_magnitude();
    return math_exp(-distsq / Superclass::m_KernelWidthSquared);
  }

  /// Evaluates \f$ K(x_{row_x},y_{row_y}) \f$.
  virtual ScalarType EvaluateKernel(const MatrixType &X, const MatrixType &Y, size_t row_x, size_t row_y) {
    ScalarType dist_squared = 0.0;
    ScalarType diff;
    for (size_t j = 0; j < X.cols(); j++) {
      diff = X(row_x, j) - Y(row_y, j);
      dist_squared += diff * diff;
    }
    return math_exp(-dist_squared / Superclass::m_KernelWidthSquared);
  }

  virtual VectorType EvaluateKernelGradient(const VectorType &x, const VectorType &y) {
    // VectorType dvec = x - y;
    ScalarType distsq = (x - y).squared_magnitude();
    ScalarType k = math_exp(-distsq / Superclass::m_KernelWidthSquared);
    // dvec *= -2.0f * k / Superclass::m_KernelWidthSquared;
    // return dvec;
    return ((x - y) * (-2.0f * k / Superclass::m_KernelWidthSquared));
  }

  /// Evaluates the gradient of \f$ K(x_{row_x},y_{row_y}) \f$ at \e \f$ x_{row_x} \f$.
  virtual VectorType EvaluateKernelGradient(const MatrixType &X, const MatrixType &Y, size_t row_x, size_t row_y) {
    VectorType result(PointDim);
    ScalarType dist_squared = 0.0;
    ScalarType diff;
    for (size_t j = 0; j < X.cols(); j++) {
      diff = X(row_x, j) - Y(row_y, j);
      result(j) = diff;
      dist_squared += diff * diff;
    }

    return (result
        * (-2.0f * math_exp(-dist_squared / Superclass::m_KernelWidthSquared) / Superclass::m_KernelWidthSquared));
  }

  // Hessian of kernel(x-y) at x
  virtual MatrixType EvaluateKernelHessian(const VectorType &x, const VectorType &y) {
    int dims = x.size();

    MatrixType xy(dims, 1, 0.0);
    ScalarType distsq = 0;
    for (int d = 0; d < dims; d++) {
      ScalarType t = x[d] - y[d];
      xy(d, 0) = t;
      distsq += t * t;
    }

    MatrixType H = xy * xy.transpose();

    ScalarType k = math_exp(-distsq / Superclass::m_KernelWidthSquared);
    H *= 4.0 * k / (Superclass::m_KernelWidthSquared * Superclass::m_KernelWidthSquared);

    for (int i = 0; i < dims; i++)
      H(i, i) -= 2.0f * k / Superclass::m_KernelWidthSquared;

    return H;
  }

  /// Evaluates the hessian of \f$ K(x_{row_x},y_{row_y}) \f$ at \e \f$ x_{row_x} \f$.
  virtual MatrixType EvaluateKernelHessian(const MatrixType &X, const MatrixType &Y, size_t row_x, size_t row_y) {
    MatrixType H(PointDim, PointDim, 0.0);
    VectorType x_minus_y(PointDim, 0.0);

    ScalarType dist_squared = 0.0;
    ScalarType diff;
    for (size_t j = 0; j < X.cols(); j++) {
      diff = X(row_x, j) - Y(row_y, j);
      x_minus_y(j) = diff;
      dist_squared += diff * diff;
    }

    for (int j = 0; j < PointDim; j++)
      for (int i = 0; i < PointDim; i++)
        H(i, j) =
            4.0 * x_minus_y(i) * x_minus_y(j) / (Superclass::m_KernelWidthSquared * Superclass::m_KernelWidthSquared);

    for (int i = 0; i < PointDim; i++)
      H(i, i) -= 2.0f / Superclass::m_KernelWidthSquared;

    ScalarType k_ij = math_exp(-dist_squared / Superclass::m_KernelWidthSquared);
    return k_ij
        * H; //(-2.0f * exp(-dist_squared / Superclass::m_KernelWidthSquared) / Superclass::m_KernelWidthSquared ) * result;
  }

  virtual MatrixType ConvolveImageFast(const MatrixType &X,const ImageTypePointer image);
  virtual std::vector<MatrixType> ConvolveGradientImageFast(const MatrixType &X, const ImageTypePointer image);
  virtual MatrixType ConvolveGradientImageFast(const MatrixType &X, const MatrixType &alpha, const ImageTypePointer img);
  virtual MatrixType Convolve(const MatrixType &X);

  /// Convolve inverse + utilities.
  MatrixType ComputeKernelMatrix(const MatrixType &Y);
  MatrixType ComputeInverseMatrix(const MatrixType &M);
  MatrixType ConvolveInverse(const MatrixType &X);

  virtual std::vector<MatrixType> ConvolveGradient(const MatrixType &X);
  virtual MatrixType ConvolveGradient(const MatrixType &X, const MatrixType &alpha);

  // Derivative of convolved weight w at direction dp
  virtual VectorType ConvolveGradient(const MatrixType &X, unsigned int k, unsigned int dp);
  virtual MatrixType ConvolveGradient(const MatrixType &X, unsigned int dim);

  virtual std::vector<std::vector<MatrixType> > ConvolveHessian(const MatrixType &X);

  // Second derivative of weight w at directions dp and dq
  virtual VectorType ConvolveHessian(const MatrixType &X, unsigned int k, unsigned int dp, unsigned int dq);
  virtual MatrixType ConvolveHessian(const MatrixType &X, unsigned int row, unsigned int col);

 protected:

}; /* class ExactKernel */

