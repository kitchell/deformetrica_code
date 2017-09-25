/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "MultiScalarInverseWishartDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
MultiScalarInverseWishartDistribution<ScalarType>
::MultiScalarInverseWishartDistribution() {
  m_DegreesOfFreedom.set_size(1);
  m_DegreesOfFreedom(0) = -1.0;

  m_ScaleVector.set_size(1);
  m_ScaleVector(0) = -1.0;
}

template<class ScalarType>
MultiScalarInverseWishartDistribution<ScalarType>
::~MultiScalarInverseWishartDistribution() {}

template<class ScalarType>
MultiScalarInverseWishartDistribution<ScalarType>
::MultiScalarInverseWishartDistribution(const MultiScalarInverseWishartDistribution &other) {
  m_DegreesOfFreedom = other.m_DegreesOfFreedom;
  m_ScaleVector = other.m_ScaleVector;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
ScalarType
MultiScalarInverseWishartDistribution<ScalarType>
::ComputeLogLikelihood(LinearVariableType const &obs) const // Input : variances.
{
  const VectorType aux = obs.vectorize();

  ScalarType out = 0.0;
  for (unsigned int k = 0; k < aux.size(); ++k) {
    out -= 0.5 * m_DegreesOfFreedom(k) * (m_ScaleVector(k) / aux(k) + std::log(aux(k)));
  }

  return out;
}

template class MultiScalarInverseWishartDistribution<double>;
