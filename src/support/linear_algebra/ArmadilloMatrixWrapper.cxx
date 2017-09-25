/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

/// Support files.
#include "ArmadilloMatrixWrapper.h"

#ifndef DEFORMETRICA_CONFIG
#include "DeformetricaConfig.h"
#endif

#ifdef USE_DOUBLE_PRECISION
#define ScalarBaseType double
#define ScalarPrecisionType float
#else
#define ScalarBaseType float
#define ScalarPrecisionType double
#endif

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator+(const ArmadilloMatrixWrapper<ScalarType> &left,
                                             const ArmadilloMatrixWrapper<ScalarType> &right) {
  return ArmadilloMatrixWrapper<ScalarType>(left.toArmadillo() + right.toArmadillo());
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator-(const ArmadilloMatrixWrapper<ScalarType> &left,
                                             const ArmadilloMatrixWrapper<ScalarType> &right) {
  return ArmadilloMatrixWrapper<ScalarType>(left.toArmadillo() - right.toArmadillo());
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> &left,
                                             const ArmadilloMatrixWrapper<ScalarType> &right) {
  return ArmadilloMatrixWrapper<ScalarType>(left.toArmadillo() * right.toArmadillo());
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarType const &rightScalar) {
  return ArmadilloMatrixWrapper<ScalarType>(leftMatrix.toArmadillo() * rightScalar);
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarPrecisionType const &rightScalar) {
  return leftMatrix * ((ScalarBaseType) rightScalar);
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(ScalarType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix) {
  return ArmadilloMatrixWrapper<ScalarType>(leftScalar * rightMatrix.toArmadillo());
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator*(ScalarPrecisionType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix) {
  return ((ScalarBaseType) leftScalar) * rightMatrix;
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(ArmadilloMatrixWrapper<ScalarType> const &leftMatrix,
                                             ArmadilloVectorWrapper<ScalarType> const &rightVector) {
  return ArmadilloVectorWrapper<ScalarType>(leftMatrix.toArmadillo() * rightVector.toArmadillo());
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarType const &rightScalar) {
  return ArmadilloMatrixWrapper<ScalarType>(leftMatrix.toArmadillo() / rightScalar);
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(const ArmadilloMatrixWrapper<ScalarType> &leftMatrix,
                                             ScalarPrecisionType const &rightScalar) {
  return leftMatrix / ((ScalarBaseType) rightScalar);
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(ScalarType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix) {
  return ArmadilloMatrixWrapper<ScalarType>(leftScalar / rightMatrix.toArmadillo());
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> operator/(ScalarPrecisionType const &leftScalar,
                                             const ArmadilloMatrixWrapper<ScalarType> &rightMatrix) {
  return ((ScalarBaseType) leftScalar) / rightMatrix;
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(unsigned N, ScalarType const &value) {
  typename ArmadilloMatrixWrapper<ScalarType>::ArmadilloMatrixType result;
  return ArmadilloMatrixWrapper<ScalarType>(result.eye(N, N) * value);
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(unsigned N, ScalarPrecisionType const &value) {
  typename ArmadilloMatrixWrapper<ScalarType>::ArmadilloMatrixType result;
  return ArmadilloMatrixWrapper<ScalarType>(result.eye(N, N) * (ScalarBaseType) value);
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> diagonal_matrix(ArmadilloVectorWrapper<ScalarType> const &values) {
  return ArmadilloMatrixWrapper<ScalarType>(arma::diagmat(values.toArmadillo()));
}
template<class ScalarType>
ScalarType trace(ArmadilloMatrixWrapper<ScalarType> const &M) {
  return arma::trace(M.toArmadillo());
}
template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> chol(ArmadilloMatrixWrapper<ScalarType> const &M) {
  return ArmadilloMatrixWrapper<ScalarType>(arma::chol(M.toArmadillo()));
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> inverse(ArmadilloMatrixWrapper<ScalarType> const &M) {
  return ArmadilloMatrixWrapper<ScalarType>(arma::inv(M.toArmadillo()));
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> inverse_sympd(ArmadilloMatrixWrapper<ScalarType> const &M) {
  // assert(arma::rcond(M.toArmadillo()) > 1e-30);
  return ArmadilloMatrixWrapper<ScalarType>(arma::inv_sympd(M.toArmadillo()));
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> eigenvalues_sym(ArmadilloMatrixWrapper<ScalarType> const &M) {
  //	arma::Col<ScalarType> eigen_values;
  typename ArmadilloVectorWrapper<ScalarType>::ArmadilloVectorType eigen_values;
  arma::eig_sym(eigen_values, M.toArmadillo());
  return ArmadilloVectorWrapper<ScalarType>(eigen_values);
}

template<class ScalarType>
std::ostream &operator<<(std::ostream &os, ArmadilloMatrixWrapper<ScalarType> const &rhs) {
  return (os << rhs.toArmadillo());
}

template<class ScalarType>
ScalarType dot_product(ArmadilloMatrixWrapper<ScalarType> const &left,
                       ArmadilloMatrixWrapper<ScalarType> const &right) {
  ScalarType result = 0;
  for (unsigned int j = 0; j < left.cols(); j++)
    result += dot_product(left.get_column(j), right.get_column(j));
  return result;
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> solve(ArmadilloMatrixWrapper<ScalarType> const &A,
                                         ArmadilloVectorWrapper<ScalarType> const &b) {
  return ArmadilloVectorWrapper<ScalarType>(arma::solve(A.toArmadillo(), b.toArmadillo()));
}

template<class ScalarType>
ArmadilloMatrixWrapper<ScalarType> solve(ArmadilloMatrixWrapper<ScalarType> const &A,
                                         ArmadilloMatrixWrapper<ScalarType> const &B) {
  return ArmadilloMatrixWrapper<ScalarType>(arma::solve(A.toArmadillo(), B.toArmadillo()));
}

/// Determinant of the input matrix \e M.
template<class ScalarType>
inline ScalarType det(ArmadilloMatrixWrapper<ScalarType> const &M) { return arma::det(M.toArmadillo()); }

/// Log determinant of the input matrix \e M.
/// \warning The determinant of \e M is assumed positive.
template<class ScalarType>
inline ScalarType log_det(ArmadilloMatrixWrapper<ScalarType> const &M) {
  ScalarType val, sign;
  arma::log_det(val, sign, M.toArmadillo());
  return val;
}

template ArmadilloMatrixWrapper<double> operator+(const ArmadilloMatrixWrapper<double> &left, const ArmadilloMatrixWrapper<double> &right) ;
template ArmadilloMatrixWrapper<double> operator-(const ArmadilloMatrixWrapper<double> &left, const ArmadilloMatrixWrapper<double> &right) ;
template ArmadilloVectorWrapper<double> operator*(ArmadilloMatrixWrapper<double> const &leftMatrix, ArmadilloVectorWrapper<double> const &rightVector);
template ArmadilloMatrixWrapper<double> operator*(const ArmadilloMatrixWrapper<double> &left, const ArmadilloMatrixWrapper<double> &right) ;
template ArmadilloMatrixWrapper<double> operator*(const ArmadilloMatrixWrapper<double> &leftMatrix, double const &rightScalar) ;
template ArmadilloMatrixWrapper<double> operator*(const ArmadilloMatrixWrapper<double> &leftMatrix, ScalarPrecisionType const &rightScalar) ;
template ArmadilloMatrixWrapper<double> operator*(double const &leftScalar, const ArmadilloMatrixWrapper<double> &rightMatrix);
template ArmadilloMatrixWrapper<double> operator*(ScalarPrecisionType const &leftScalar, const ArmadilloMatrixWrapper<double> &rightMatrix);
template ArmadilloMatrixWrapper<double> operator/(const ArmadilloMatrixWrapper<double> &leftMatrix, double const &rightScalar) ;
template ArmadilloMatrixWrapper<double> operator/(const ArmadilloMatrixWrapper<double> &leftMatrix, ScalarPrecisionType const &rightScalar) ;
template ArmadilloMatrixWrapper<double> operator/(double const &leftScalar, const ArmadilloMatrixWrapper<double> &rightMatrix);
template ArmadilloMatrixWrapper<double> operator/(ScalarPrecisionType const &leftScalar, const ArmadilloMatrixWrapper<double> &rightMatrix);
template double trace<double>(ArmadilloMatrixWrapper<double> const &M);
template ArmadilloVectorWrapper<double> solve<double>(ArmadilloMatrixWrapper<double> const &A, ArmadilloVectorWrapper<double> const &b);
template ArmadilloMatrixWrapper<double> solve<double>(ArmadilloMatrixWrapper<double> const &A, ArmadilloMatrixWrapper<double> const &B);
template ArmadilloMatrixWrapper<double> diagonal_matrix<double>(unsigned N, double const &value);
template ArmadilloMatrixWrapper<double> diagonal_matrix<double>(unsigned N, ScalarPrecisionType const &value);
template ArmadilloMatrixWrapper<double> diagonal_matrix<double>(ArmadilloVectorWrapper<double> const &values);
template ArmadilloMatrixWrapper<double> chol<double>(ArmadilloMatrixWrapper<double> const &M);
template ArmadilloMatrixWrapper<double> inverse<double>(ArmadilloMatrixWrapper<double> const &M);
template ArmadilloMatrixWrapper<double> inverse_sympd<double>(ArmadilloMatrixWrapper<double> const &M);
template ArmadilloVectorWrapper<double> eigenvalues_sym<double>(ArmadilloMatrixWrapper<double> const &M);
template std::ostream &operator<<(std::ostream &os, ArmadilloMatrixWrapper<double> const &rhs);
template double dot_product<double>(ArmadilloMatrixWrapper<double> const &left, ArmadilloMatrixWrapper<double> const &right);
template double det<double>(ArmadilloMatrixWrapper<double> const &M);
template double log_det<double>(ArmadilloMatrixWrapper<double> const &M);

#undef ScalarBaseType
#undef ScalarPrecisionType
