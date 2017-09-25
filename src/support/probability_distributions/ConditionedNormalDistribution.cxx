/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "ConditionedNormalDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
ConditionedNormalDistribution<ScalarType>
::ConditionedNormalDistribution() : Superclass() {}

template<class ScalarType>
ConditionedNormalDistribution<ScalarType>
::~ConditionedNormalDistribution() {}

template<class ScalarType>
ConditionedNormalDistribution<ScalarType>
::ConditionedNormalDistribution(const ConditionedNormalDistribution &other) : Superclass(other) {}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
VectorType
ConditionedNormalDistribution<ScalarType>
::Sample() const {
  VectorType out;
  unsigned int n = 0;

  for (; n < 100; ++n) {
    /// Sampling according to a simple normal distribution.
    out = Superclass::Sample();

    /// Clamps the sampled values between 0 and 1.
    for (unsigned int k = 0; k < out.size(); ++k) {
      if (out(k) < 0.0) { out(k) = 0.0; }
      if (out(k) > 1.0) { out(k) = 1.0; }
    }

    /// Normalization of the vector, which must sum to 1.
    ScalarType sum = out.sum();
    if (sum > 1e-10) { return out / sum; }
  }

  std::cerr << "Warning : in ConditionedNormalDistribution::Sample(), after 100 tries, no suitable output "
      "has been found. Defaults to the mean of the normal distribution." << std::endl;
  return Superclass::m_Mean;
}

template class ConditionedNormalDistribution<double>;
