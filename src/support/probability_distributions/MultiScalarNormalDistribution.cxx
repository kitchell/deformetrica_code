/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "MultiScalarNormalDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
MultiScalarNormalDistribution<ScalarType>
::MultiScalarNormalDistribution() {
  m_Variance = -1.0;
  m_VarianceSqrt = -1.0;
  m_VarianceInverse = -1.0;
  m_VarianceLog = -1.0;
}

template<class ScalarType>
MultiScalarNormalDistribution<ScalarType>
::~MultiScalarNormalDistribution() {}

template<class ScalarType>
MultiScalarNormalDistribution<ScalarType>
::MultiScalarNormalDistribution(const MultiScalarNormalDistribution &other) : Superclass(other) {
  m_Variance = other.m_Variance;
  m_VarianceSqrt = other.m_VarianceSqrt;
  m_VarianceInverse = other.m_VarianceInverse;
  m_VarianceLog = other.m_VarianceLog;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
VectorType
MultiScalarNormalDistribution<ScalarType>
::Sample() const {
  std::random_device rd;
  std::mt19937_64 generator(rd());
  std::normal_distribution<ScalarType> stdNormal(0, 1);

  const unsigned int N = Superclass::m_Mean.size();
  VectorType aux(N);

  for (unsigned int i = 0; i < N; ++i)
    aux(i) = stdNormal(generator);

  return Superclass::m_Mean + m_VarianceSqrt * aux;
}

template<class ScalarType>
ScalarType
MultiScalarNormalDistribution<ScalarType>
::ComputeLogLikelihood(VectorType const &obs,
                       ScalarType const &temperature) const {
  return -0.5 * (m_VarianceInverse / temperature) * (Superclass::m_Mean - obs).sum_of_squares()
      - 0.5 * Superclass::m_Mean.size() * (m_VarianceLog + std::log(temperature));
}

template class MultiScalarNormalDistribution<double>;
