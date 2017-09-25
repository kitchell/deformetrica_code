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
#include "ArmadilloVectorWrapper.h"

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
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> &left,
                                             const ArmadilloVectorWrapper<ScalarType> &right) {
  return ArmadilloVectorWrapper<ScalarType>(left.toArmadillo() - right.toArmadillo());
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> &left,
                                             const ArmadilloVectorWrapper<ScalarType> &right) {
  return ArmadilloVectorWrapper<ScalarType>(left.toArmadillo() + right.toArmadillo());
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator%(const ArmadilloVectorWrapper<ScalarType> &left,
                                             const ArmadilloVectorWrapper<ScalarType> &right) {
  return ArmadilloVectorWrapper<ScalarType>(left.toArmadillo() % right.toArmadillo());
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarType &rightScalar) {
  return ArmadilloVectorWrapper<ScalarType>(leftVector.toArmadillo() + rightScalar);
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator+(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarPrecisionType &rightScalar) {
  return leftVector + ((ScalarBaseType) rightScalar);
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarType &rightScalar) {
  return ArmadilloVectorWrapper<ScalarType>(leftVector.toArmadillo() - rightScalar);
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator-(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarPrecisionType &rightScalar) {
  return leftVector - ((ScalarBaseType) rightScalar);
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarType &rightScalar) {
  return ArmadilloVectorWrapper<ScalarType>(leftVector.toArmadillo() / rightScalar);
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ArmadilloVectorWrapper<ScalarType> &leftVector,
                                             const ScalarPrecisionType &rightScalar) {
  return leftVector / ((ScalarBaseType) rightScalar);
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ScalarType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector) {
  unsigned int size = rightVector.size();
  ArmadilloVectorWrapper<ScalarType> out(size);
  for (unsigned int i = 0; i < size; ++i) { out[i] = leftScalar / rightVector[i]; }
  return out;
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator/(const ScalarPrecisionType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector) {
  return ((ScalarBaseType) leftScalar) / rightVector;
}

template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(const ScalarType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector) {
  return ArmadilloVectorWrapper<ScalarType>(leftScalar * rightVector.toArmadillo());
}
template<class ScalarType>
ArmadilloVectorWrapper<ScalarType> operator*(const ScalarPrecisionType &leftScalar,
                                             const ArmadilloVectorWrapper<ScalarType> &rightVector) {
  return ((ScalarBaseType) leftScalar) * rightVector;
}

template<class ScalarType>
std::ostream &operator<<(std::ostream &os, ArmadilloVectorWrapper<ScalarType> const &rhs) {
  return (os << rhs.toArmadillo().t());
}

template<class ScalarType>
inline ScalarType dot_product(ArmadilloVectorWrapper<ScalarType> const &v1,
                              ArmadilloVectorWrapper<ScalarType> const &v2) {
  return arma::dot(v1.toArmadillo(), v2.toArmadillo());
}

template<class ScalarType>
inline ArmadilloVectorWrapper<ScalarType> cross_3d(ArmadilloVectorWrapper<ScalarType> const &v1,
                                                   ArmadilloVectorWrapper<ScalarType> const &v2) {
  return ArmadilloVectorWrapper<ScalarType>(arma::cross(v1.toArmadillo(), v2.toArmadillo()));
}


template class ArmadilloVectorWrapper<float>;
template class ArmadilloVectorWrapper<double>;
template ArmadilloVectorWrapper<double> operator%(const ArmadilloVectorWrapper<double> &left, const ArmadilloVectorWrapper<double> &right) ;
template ArmadilloVectorWrapper<double> operator+(const ArmadilloVectorWrapper<double> &left, const ArmadilloVectorWrapper<double> &right) ;
template ArmadilloVectorWrapper<double> operator+(const ArmadilloVectorWrapper<double> &leftVector, const double &rightScalar) ;
template ArmadilloVectorWrapper<double> operator+(const ArmadilloVectorWrapper<double> &leftVector, const ScalarPrecisionType &rightScalar) ;
template ArmadilloVectorWrapper<double> operator-(const ArmadilloVectorWrapper<double> &left, const ArmadilloVectorWrapper<double> &right) ;
template ArmadilloVectorWrapper<double> operator-(const ArmadilloVectorWrapper<double> &leftVector, const double &rightScalar) ;
template ArmadilloVectorWrapper<double> operator-(const ArmadilloVectorWrapper<double> &leftVector, const ScalarPrecisionType &rightScalar) ;
template ArmadilloVectorWrapper<double> operator/(const ArmadilloVectorWrapper<double> &leftVector, const double &rightScalar) ;
template ArmadilloVectorWrapper<double> operator/(const ArmadilloVectorWrapper<double> &leftVector, const ScalarPrecisionType &rightScalar) ;
template ArmadilloVectorWrapper<double> operator/(const double &leftScalar, const ArmadilloVectorWrapper<double> &rightVector);
template ArmadilloVectorWrapper<double> operator/(const ScalarPrecisionType &leftScalar, const ArmadilloVectorWrapper<double> &rightVector) ;
template ArmadilloVectorWrapper<double> operator*(const double &leftScalar, const ArmadilloVectorWrapper<double> &rightVector) ;
template ArmadilloVectorWrapper<double> operator*(const ScalarPrecisionType &leftScalar, const ArmadilloVectorWrapper<double> &rightVector) ;
template std::ostream &operator<<(std::ostream &os, ArmadilloVectorWrapper<double> const &rhs) ;
template double dot_product(ArmadilloVectorWrapper<double> const &v1, ArmadilloVectorWrapper<double> const &v2) ;
template ArmadilloVectorWrapper<double> cross_3d(ArmadilloVectorWrapper<double> const &v1, ArmadilloVectorWrapper<double> const &v2) ;

#undef ScalarBaseType
#undef ScalarPrecisionType