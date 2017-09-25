/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "NormalDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
NormalDistribution<ScalarType>
::NormalDistribution() {
  m_Covariance.set_size(1, 1);
  m_Covariance(0, 0) = 1.0;
  m_CovarianceSqrt.set_size(1, 1);
  m_CovarianceSqrt(0, 0) = 1.0;
  m_CovarianceInverse.set_size(1, 1);
  m_CovarianceInverse(0, 0) = 1.0;
}

template<class ScalarType>
NormalDistribution<ScalarType>
::~NormalDistribution() {}

template<class ScalarType>
NormalDistribution<ScalarType>
::NormalDistribution(const NormalDistribution &other) : Superclass(other) {
  m_Covariance = other.m_Covariance;
  m_CovarianceSqrt = other.m_CovarianceSqrt;
  m_CovarianceInverse = other.m_CovarianceInverse;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
VectorType
NormalDistribution<ScalarType>
::Sample() const {
  std::random_device rd;
  std::mt19937_64 generator(rd());
  std::normal_distribution<ScalarType> stdNormal(0, 1);

  const unsigned int N = Superclass::m_Mean.size();
  VectorType aux(N);

  for (int i = 0; i < N; i++)
    aux(i) = stdNormal(generator);

  return Superclass::m_Mean + m_CovarianceSqrt * aux;
}

template<class ScalarType>
ScalarType
NormalDistribution<ScalarType>
::ComputeLogLikelihood(VectorType const &obs,
                       ScalarType const &temperature) const {
  return -0.5 * dot_product(obs - Superclass::m_Mean, (m_CovarianceInverse / temperature) * (obs - Superclass::m_Mean))
      - 0.5 * (m_CovarianceLogDeterminant + Superclass::m_Mean.size() * std::log(temperature));
}

template class NormalDistribution<double>;
