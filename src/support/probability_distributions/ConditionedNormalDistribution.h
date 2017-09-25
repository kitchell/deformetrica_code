/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _ConditionedNormalDistribution_h
#define _ConditionedNormalDistribution_h

/// Class file.
#include "NormalDistribution.h"

/**
 *  \brief      ConditionedNormalDistribution class
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    Allows to draw samples from a normal distribution, determined by its mean and variance.
 */
template<class ScalarType>
class ConditionedNormalDistribution : public NormalDistribution<ScalarType> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract probability distribution type.
  typedef NormalDistribution<ScalarType> Superclass;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructors.
  ConditionedNormalDistribution(); // Standard normal distribution.

  /// Copy constructor.
  ConditionedNormalDistribution(const ConditionedNormalDistribution &other);

  /// Makes a copy of the object.
  virtual ConditionedNormalDistribution *Clone() { return new ConditionedNormalDistribution(*this); }

  /// Destructor.
  virtual ~ConditionedNormalDistribution();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Samples from the normal distribution to generate a realization.
  virtual VectorType Sample() const;

}; /* class ConditionedNormalDistribution */


//#ifndef MU_MANUAL_INSTANTIATION
//#include "ConditionedNormalDistribution.txx"
//#endif

#endif /* _ConditionedNormalDistribution_h */
