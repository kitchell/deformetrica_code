/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

//#ifndef _MatrixListWrapper_friend_h
//#define _MatrixListWrapper_friend_h
//
//#ifndef _MatrixListWrapper_h
//#error Do not include MatrixListWrapper_friend.h : include MatrixListWrapper.h instead
//#endif

#include "MatrixListWrapper.h"

//
// Left matrix list / Right matrix list :
//
template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator+(const MatrixListWrapper<ScalarType> &left,
                                               const MatrixListWrapper<ScalarType> &right) {
  std::vector<ArmadilloMatrixWrapper<ScalarType>> rawMatrixListType(left.m_RawMatrixList.size());
  for (unsigned int k = 0; k < left.m_RawMatrixList.size(); k++)
    rawMatrixListType[k] = left.m_RawMatrixList[k] + right.m_RawMatrixList[k];
  return MatrixListWrapper<ScalarType>(rawMatrixListType);
}

template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator-(const MatrixListWrapper<ScalarType> &left,
                                               const MatrixListWrapper<ScalarType> &right) {
  std::vector<ArmadilloMatrixWrapper<ScalarType>> rawMatrixListType(left.m_RawMatrixList.size());
  for (unsigned int k = 0; k < left.m_RawMatrixList.size(); k++)
    rawMatrixListType[k] = left.m_RawMatrixList[k] - right.m_RawMatrixList[k];
  return MatrixListWrapper<ScalarType>(rawMatrixListType);
}

//
// Left-Right matrix list / scalar :
//
template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator*(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                               ScalarType const &rightScalar) {
  std::vector<ArmadilloMatrixWrapper<ScalarType>> rawMatrixListType(leftMatrixListType.m_RawMatrixList.size());
  for (unsigned int k = 0; k < leftMatrixListType.m_RawMatrixList.size(); k++)
    rawMatrixListType[k] = leftMatrixListType.m_RawMatrixList[k] * rightScalar;
  return MatrixListWrapper<ScalarType>(rawMatrixListType);
}

template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator*(ScalarType const &leftScalar,
                                               const MatrixListWrapper<ScalarType> &rightMatrixListType) {
  std::vector<ArmadilloMatrixWrapper<ScalarType>> rawMatrixListType(rightMatrixListType.m_RawMatrixList.size());
  for (unsigned int k = 0; k < rightMatrixListType.m_RawMatrixList.size(); k++)
    rawMatrixListType[k] = leftScalar * rightMatrixListType.m_RawMatrixList[k];
  return MatrixListWrapper<ScalarType>(rawMatrixListType);
}

template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator/(const MatrixListWrapper<ScalarType> &leftMatrixListType,
                                               ScalarType const &rightScalar) {
  std::vector<ArmadilloMatrixWrapper<ScalarType>> rawMatrixListType(leftMatrixListType.m_RawMatrixList.size());
  for (unsigned int k = 0; k < leftMatrixListType.m_RawMatrixList.size(); k++)
    rawMatrixListType[k] = leftMatrixListType.m_RawMatrixList[k] / rightScalar;
  return MatrixListWrapper<ScalarType>(rawMatrixListType);
}

//
// Left-Right matrix list / int :
//
template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator*(const MatrixListWrapper<ScalarType> &left,
                                               int const &right) {
  return operator*(left, ScalarType(right));
}

template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator*(int const &left,
                                               const MatrixListWrapper<ScalarType> &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline MatrixListWrapper<ScalarType> operator/(const MatrixListWrapper<ScalarType> &left,
                                               int const &right) {
  return operator*(left, 1 / right);
}

template
MatrixListWrapper<double> operator+(const MatrixListWrapper<double> &left,
                                    const MatrixListWrapper<double> &right);
template
MatrixListWrapper<double> operator-(const MatrixListWrapper<double> &left,
                                    const MatrixListWrapper<double> &right);

template
MatrixListWrapper<double> operator*(const MatrixListWrapper<double> &leftMatrixListType,
                                    double const &rightScalar);
template
MatrixListWrapper<double> operator*(double const &leftScalar,
                                    const MatrixListWrapper<double> &rightMatrixListType);
template
MatrixListWrapper<double> operator/(const MatrixListWrapper<double> &leftMatrixListType,
                                    double const &rightScalar);

template
MatrixListWrapper<double> operator*(const MatrixListWrapper<double> &leftMatrixListType,
                                    int const &rightScalar);
template
MatrixListWrapper<double> operator*(int const &leftScalar,
                                    const MatrixListWrapper<double> &rightMatrixListType);
template
MatrixListWrapper<double> operator/(const MatrixListWrapper<double> &leftMatrixListType,
                                    int const &rightScalar);

