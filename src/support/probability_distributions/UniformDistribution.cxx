/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "UniformDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
UniformDistribution<ScalarType>
::UniformDistribution() {
  m_LowerBounds.set_size(1);
  m_LowerBounds(0) = 0.0;
  m_UpperBounds.set_size(1);
  m_UpperBounds(0) = 1.0;
}

template<class ScalarType>
UniformDistribution<ScalarType>
::UniformDistribution(VectorType const &lowers,
                      VectorType const &uppers) : m_LowerBounds(lowers), m_UpperBounds(uppers) {}

template<class ScalarType>
UniformDistribution<ScalarType>
::~UniformDistribution() {}

template<class ScalarType>
UniformDistribution<ScalarType>
::UniformDistribution(const UniformDistribution &other) {
  m_LowerBounds = other.m_LowerBounds;
  m_UpperBounds = other.m_UpperBounds;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
VectorType
UniformDistribution<ScalarType>
::Sample() const {
  std::random_device rd;
  std::mt19937_64 generator(rd());

  const unsigned int N = m_LowerBounds.size();
  VectorType out(N);

  for (unsigned int k = 0; k < m_LowerBounds.size(); ++k) {
    std::uniform_real_distribution<ScalarType> distribution(m_LowerBounds(k), m_UpperBounds(k));
    out(k) = distribution(generator);
  }

  return out;
}

template<class ScalarType>
ScalarType
UniformDistribution<ScalarType>
::ComputeLogLikelihood(LinearVariableType const &obs) const {
  const unsigned int N = m_LowerBounds.size();
  assert(N == obs.n_elem());

  VectorType aux = obs.vectorize();

  ScalarType out = 0.0;
  for (unsigned int k = 0; k < m_LowerBounds.size(); ++k) {
    if (aux(k) >= m_LowerBounds(k) && aux(k) <= m_UpperBounds(k))
      out -= std::log(m_UpperBounds(k) - m_LowerBounds(k));
  }

  return out;
}

template class UniformDistribution<double>;