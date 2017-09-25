/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "LinearVariablesMapWrapper.h"

//
// Left linear variable map / Right linear variable map :
//
template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator+(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       LinearVariablesMapWrapper<ScalarType> const &right) {
  LinearVariablesMapWrapper<ScalarType> out = left;
  return out += right;
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator-(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       LinearVariablesMapWrapper<ScalarType> const &right) {
  LinearVariablesMapWrapper<ScalarType> out = left;
  return out -= right;
}

//
// Left linear variable map / Right scalar :
//
template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       ScalarType const &right) {
  LinearVariablesMapWrapper<ScalarType> out = left;
  return out *= right;
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(ScalarType const &left,
                                                       LinearVariablesMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator/(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       ScalarType const &right) {
  return operator*(left, 1 / right);
}

//
// Left linear variable map / Right int :
//
template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       int const &right) {
  LinearVariablesMapWrapper<ScalarType> out = left;
  return out *= right;
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(int const &left,
                                                       LinearVariablesMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator/(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       int const &right) {
  return operator*(left, 1 / right);
}

//
// Left linear variable map / Right unsigned int :
//
template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       unsigned int const &right) {
  LinearVariablesMapWrapper<ScalarType> out = left;
  return out *= right;
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(unsigned int const &left,
                                                       LinearVariablesMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator/(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       unsigned int const &right) {
  return operator*(left, 1 / right);
}

//
// Left linear variable map / Right vector :
//
template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(LinearVariablesMapWrapper<ScalarType> const &left,
                                                       ArmadilloVectorWrapper<ScalarType> const &right) {
  assert(left.size() == right.size());
  LinearVariablesMapWrapper<ScalarType> out = left;
  unsigned int ind = 0;
  for (auto it = out.m_RawLinearVariablesMap.begin(); it != out.m_RawLinearVariablesMap.end(); ++it, ++ind) {
    for (unsigned int k = 0; k < it->second.size(); ++k) {
      it->second[k] *= right[ind];
    }
  }
  return out;
}

template<class ScalarType>
inline LinearVariablesMapWrapper<ScalarType> operator*(ArmadilloVectorWrapper<ScalarType> const &left,
                                                       LinearVariablesMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}


//*****************Define TEMPLATE

//
// Left linear variable map / Right linear variable map :
//
template
LinearVariablesMapWrapper<double> operator+(LinearVariablesMapWrapper<double> const &left,
                                            LinearVariablesMapWrapper<double> const &right);
template
LinearVariablesMapWrapper<double> operator-(LinearVariablesMapWrapper<double> const &left,
                                            LinearVariablesMapWrapper<double> const &right);

template
LinearVariablesMapWrapper<double> operator*(LinearVariablesMapWrapper<double> const &left,
                                            ArmadilloVectorWrapper<double> const &right);
template
LinearVariablesMapWrapper<double> operator*(ArmadilloVectorWrapper<double> const &left,
                                            LinearVariablesMapWrapper<double> const &right);

template
LinearVariablesMapWrapper<double> operator*(LinearVariablesMapWrapper<double> const &left, double const &right);
template
LinearVariablesMapWrapper<double> operator*(double const &left, LinearVariablesMapWrapper<double> const &right);
template
LinearVariablesMapWrapper<double> operator/(LinearVariablesMapWrapper<double> const &left, double const &right);

template
LinearVariablesMapWrapper<double> operator*(LinearVariablesMapWrapper<double> const &left, int const &right);
template
LinearVariablesMapWrapper<double> operator*(int const &left, LinearVariablesMapWrapper<double> const &right);
template
LinearVariablesMapWrapper<double> operator/(LinearVariablesMapWrapper<double> const &left, int const &right);

template
LinearVariablesMapWrapper<double> operator*(LinearVariablesMapWrapper<double> const &left, unsigned int const &right);
template
LinearVariablesMapWrapper<double> operator*(unsigned int const &left, LinearVariablesMapWrapper<double> const &right);
template
LinearVariablesMapWrapper<double> operator/(LinearVariablesMapWrapper<double> const &left, unsigned int const &right);


