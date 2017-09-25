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

/// Class file.
#include "AbstractProbabilityDistribution.h"
#include "LinearAlgebra.h"
#include <math.h>

using namespace def::algebra;

/**
 *  \brief      AutomaticRelevanceDeterminationDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Implementation of an ARD distribution with parameter gamma (p(W) = prod_{i=1}^N (gamma_i/2pi)
 */
template<class ScalarType>
class AutomaticRelevanceDeterminationDistribution : public AbstractProbabilityDistribution<ScalarType> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract probability distribution type.
  typedef AbstractProbabilityDistribution<ScalarType> Superclass;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructors.
  AutomaticRelevanceDeterminationDistribution();

  /// Copy constructor.
  AutomaticRelevanceDeterminationDistribution(const AutomaticRelevanceDeterminationDistribution &other);

  /// Makes a copy of the object.
  virtual AutomaticRelevanceDeterminationDistribution *Clone() { return new AutomaticRelevanceDeterminationDistribution(*this); }

  /// Destructor.
  virtual ~AutomaticRelevanceDeterminationDistribution();


  /// Samples from the normal distribution to generate a realization.
  virtual VectorType Sample() const {throw(std::runtime_error("Not implemented for ARD"));};

  /// Variation featuring a temperature parameter (for the McmcSaem algorithm). Vectorizes \e obs in a first step.
  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs,
                                          ScalarType const &temperature) const {
    MatrixType aux = recast<MatrixType>(obs);
    return ComputeLogLikelihood(aux, temperature);
  }

  /// Variation featuring a temperature parameter (for the McmcSaem algorithm).
  virtual ScalarType ComputeLogLikelihood(MatrixType const &obs,
                                          ScalarType const &temperature) const;

  virtual ScalarType ComputeLogLikelihood(LinearVariableType const &obs) const {return ComputeLogLikelihood(obs, 1.);};


  void SetGamma(VectorType newGamma){m_Gamma = newGamma;};

  VectorType GetGamma() const {return m_Gamma;}

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///Gamma parameters for the ARD distribution
  VectorType m_Gamma;

}; /* class AutomaticRelevanceDeterminationDistribution */
