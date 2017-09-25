/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

/// Support files.
#include "LinearAlgebra.h"

/// Librairies.
#include <random>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace def::algebra;

/**
 *  \brief      AbstractProbabilityDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples, and to compute the likelihood of an observation.
 */
template<class ScalarType>
class AbstractProbabilityDistribution {
 public:


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  AbstractProbabilityDistribution();

  /// Copy constructor.
  AbstractProbabilityDistribution(const AbstractProbabilityDistribution &other);

  /// Makes a copy of the object.
  virtual AbstractProbabilityDistribution *Clone() = 0;

  /// Destructor.
  virtual ~AbstractProbabilityDistribution();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples from the distribution to generate a realization.
  virtual VectorType Sample() const = 0;

  /// Computes the log-likelihood of an observation in LinearVariable format.
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs) const = 0;
  /// Optional variation, featuring a temperature parameter (for the McmcSaem algorithm).
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs,
                                          ScalarType const &temperature) const {
    return ComputeLogLikelihood(obs);
  }

}; /* class AbstractProbabilityDistribution */


