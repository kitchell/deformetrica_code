/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "LinearVariableMapWrapper.h"

//
// Left linear variable map / Right linear variable map :
//
template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator+(LinearVariableMapWrapper<ScalarType> const &left,
                                                      LinearVariableMapWrapper<ScalarType> const &right) {
  typedef typename LinearVariableMapWrapper<ScalarType>::RawLinearVariableMapType RawLinearVariableMapType;
  typedef typename LinearVariableMapWrapper<ScalarType>::const_iterator const_iterator;
  addition_visitor<ScalarType> visitor;
  RawLinearVariableMapType rawList;
  for (const_iterator it = left.m_RawLinearVariableMap.begin(); it != left.m_RawLinearVariableMap.end(); ++it)
    rawList[it->first] = boost::apply_visitor(visitor, it->second.GetRawLinearVariable(),
                                              right.m_RawLinearVariableMap.at(it->first).GetRawLinearVariable());
  return LinearVariableMapWrapper<ScalarType>(rawList);
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator-(LinearVariableMapWrapper<ScalarType> const &left,
                                                      LinearVariableMapWrapper<ScalarType> const &right) {
  typedef typename LinearVariableMapWrapper<ScalarType>::RawLinearVariableMapType RawLinearVariableMapType;
  typedef typename LinearVariableMapWrapper<ScalarType>::const_iterator const_iterator;
  subtraction_visitor<ScalarType> visitor;
  RawLinearVariableMapType rawList;
  for (const_iterator it = left.m_RawLinearVariableMap.begin(); it != left.m_RawLinearVariableMap.end(); ++it)
    rawList[it->first] = boost::apply_visitor(visitor, it->second.GetRawLinearVariable(),
                                              right.m_RawLinearVariableMap.at(it->first).GetRawLinearVariable());
  return LinearVariableMapWrapper<ScalarType>(rawList);
}

//
// Left linear variable map / Right scalar :
//
template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                                      ScalarType const &right) {
  typedef typename LinearVariableMapWrapper<ScalarType>::RawLinearVariableMapType RawLinearVariableMapType;
  typedef typename LinearVariableMapWrapper<ScalarType>::const_iterator const_iterator;
  right_scalar_multiplication_visitor<ScalarType> visitor(right);
  RawLinearVariableMapType rawList;
  for (const_iterator it = left.m_RawLinearVariableMap.begin(); it != left.m_RawLinearVariableMap.end(); ++it)
    rawList[it->first] = boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
  return LinearVariableMapWrapper<ScalarType>(rawList);
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(ScalarType const &left,
                                                      LinearVariableMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator/(LinearVariableMapWrapper<ScalarType> const &left,
                                                      ScalarType const &right) {
  return operator*(left, 1 / right);
}

//
// Left linear variable map / Right int :
//
template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                                      int const &right) {
  LinearVariableMapWrapper<ScalarType> out = left;
  return out *= right;
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(int const &left,
                                                      LinearVariableMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator/(LinearVariableMapWrapper<ScalarType> const &left,
                                                      int const &right) {
  return operator*(left, 1 / right);
}

//
// Left linear variable map / Right unsigned int :
//
template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                                      unsigned int const &right) {
  LinearVariableMapWrapper<ScalarType> out = left;
  return out *= right;
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(unsigned int const &left,
                                                      LinearVariableMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator/(LinearVariableMapWrapper<ScalarType> const &left,
                                                      unsigned int const &right) {
  return operator*(left, 1 / right);
}

//
// Left linear variable map / Right vector :
//
template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(LinearVariableMapWrapper<ScalarType> const &left,
                                                      ArmadilloVectorWrapper<ScalarType> const &right) {
  assert(left.size() == right.size());
  typedef typename LinearVariableMapWrapper<ScalarType>::RawLinearVariableMapType RawLinearVariableMapType;
  typedef typename LinearVariableMapWrapper<ScalarType>::const_iterator const_iterator;
  RawLinearVariableMapType rawList;
  unsigned int k = 0;
  for (const_iterator it = left.m_RawLinearVariableMap.begin(); it != left.m_RawLinearVariableMap.end(); ++it, ++k) {
    right_scalar_multiplication_visitor<ScalarType> visitor(right(k));
    rawList[it->first] = boost::apply_visitor(visitor)(it->second.GetRawLinearVariable());
  }
  return LinearVariableMapWrapper<ScalarType>(rawList);
}

template<class ScalarType>
inline LinearVariableMapWrapper<ScalarType> operator*(ArmadilloVectorWrapper<ScalarType> const &left,
                                                      LinearVariableMapWrapper<ScalarType> const &right) {
  return operator*(right, left);
}


//*****************Define TEMPLATE

//
// Left linear variable map / Right linear variable map :
//
template
LinearVariableMapWrapper<double> operator+(LinearVariableMapWrapper<double> const &left,
                                           LinearVariableMapWrapper<double> const &right);
template
LinearVariableMapWrapper<double> operator-(LinearVariableMapWrapper<double> const &left,
                                           LinearVariableMapWrapper<double> const &right);

template
LinearVariableMapWrapper<double> operator*(LinearVariableMapWrapper<double> const &left,
                                           ArmadilloVectorWrapper<double> const &right);
template
LinearVariableMapWrapper<double> operator*(ArmadilloVectorWrapper<double> const &left,
                                           LinearVariableMapWrapper<double> const &right);

template
LinearVariableMapWrapper<double> operator*(LinearVariableMapWrapper<double> const &left, double const &right);
template
LinearVariableMapWrapper<double> operator*(double const &left, LinearVariableMapWrapper<double> const &right);
template
LinearVariableMapWrapper<double> operator/(LinearVariableMapWrapper<double> const &left, double const &right);

template
LinearVariableMapWrapper<double> operator*(LinearVariableMapWrapper<double> const &left, int const &right);
template
LinearVariableMapWrapper<double> operator*(int const &left, LinearVariableMapWrapper<double> const &right);
template
LinearVariableMapWrapper<double> operator/(LinearVariableMapWrapper<double> const &left, int const &right);

template
LinearVariableMapWrapper<double> operator*(LinearVariableMapWrapper<double> const &left, unsigned int const &right);
template
LinearVariableMapWrapper<double> operator*(unsigned int const &left, LinearVariableMapWrapper<double> const &right);
template
LinearVariableMapWrapper<double> operator/(LinearVariableMapWrapper<double> const &left, unsigned int const &right);



