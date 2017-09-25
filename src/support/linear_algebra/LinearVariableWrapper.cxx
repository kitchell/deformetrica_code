/****************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

//#ifndef _LinearVariableWrapper_friend_h
//#define _LinearVariableWrapper_friend_h
//
//#ifndef _LinearVariableWrapper_h
//#error Do not include LinearVariableWrapper_friend.h : include LinearVariableWrapper.h instead
//#endif

#include "LinearVariableWrapper.h"

//
// Left linear variable / Right linear variable :
//
template<class ScalarType>
inline LinearVariableWrapper<ScalarType> operator+(LinearVariableWrapper<ScalarType> const &left,
                                                   LinearVariableWrapper<ScalarType> const &right) {
  addition_visitor<ScalarType> visitor;
  return LinearVariableWrapper<ScalarType>
      (boost::apply_visitor(visitor, left.m_RawLinearVariable, right.m_RawLinearVariable));
}

template<class ScalarType>
inline LinearVariableWrapper<ScalarType> operator-(LinearVariableWrapper<ScalarType> const &left,
                                                   LinearVariableWrapper<ScalarType> const &right) {
  subtraction_visitor<ScalarType> visitor;
  return LinearVariableWrapper<ScalarType>
      (boost::apply_visitor(visitor)(left.m_RawLinearVariable, right.m_RawLinearVariable));
}

//
// Left linear variable / Right scalar :
//
template<class ScalarType>
inline LinearVariableWrapper<ScalarType> operator*(LinearVariableWrapper<ScalarType> const &left,
                                                   ScalarType const &right) {
  right_scalar_multiplication_visitor<ScalarType> visitor(right);
  return LinearVariableWrapper<ScalarType>(boost::apply_visitor(visitor)(left.m_RawLinearVariable));
}

template<class ScalarType>
inline LinearVariableWrapper<ScalarType> operator*(ScalarType const &left,
                                                   LinearVariableWrapper<ScalarType> const &right) {
  return operator*(right, left);
}

template<class ScalarType>
inline LinearVariableWrapper<ScalarType> operator/(LinearVariableWrapper<ScalarType> const &left,
                                                   ScalarType const &right) {
  return operator*(left, 1 / right);
}

template
LinearVariableWrapper<double> operator+(LinearVariableWrapper<double> const &left,
                                        LinearVariableWrapper<double> const &right);
template
LinearVariableWrapper<double> operator-(LinearVariableWrapper<double> const &left,
                                        LinearVariableWrapper<double> const &right);

template
LinearVariableWrapper<double> operator*(LinearVariableWrapper<double> const &left,
                                        double const &right);
template
LinearVariableWrapper<double> operator*(double const &left,
                                        LinearVariableWrapper<double> const &right);
template
LinearVariableWrapper<double> operator/(LinearVariableWrapper<double> const &left,
                                        double const &right);


//#endif /* _LinearVariableWrapper_friend_h */
