/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah. All rights reserved. This file is     *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

/// Class file.
#include "AbstractEstimator.h"

/// Core files.
#include "LinearAlgebra.h"

using namespace def::algebra;

/**
 *	\brief      PowellsMethod object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A model estimator is an algorithm which updates the fixed effects of a statistical model.
 */
template<class ScalarType, unsigned int Dimension>
class PowellsMethod : public AbstractEstimator<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Model estimator type.
  typedef AbstractEstimator<ScalarType, Dimension> Superclass;

  /// Abstract statistical model type.
  typedef typename Superclass::StatisticalModelType StatisticalModelType;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  PowellsMethod();

  /// Copy constructor.
  PowellsMethod(const PowellsMethod &other);

  /// Makes a copy of the object.
  virtual PowellsMethod *Clone() { return new PowellsMethod(*this); }

  /// Destructor.
  virtual ~PowellsMethod();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the maximum number of iterations in the line search procedure.
  unsigned int GetMaxLineSearchIterations() const { return m_MaxLineSearchIterations; }
  /// Sets the maximum number of iterations in the line search procedure to \e n.
  void SetMaxLineSearchIterations(const unsigned int n) { m_MaxLineSearchIterations = n; }

  /// Sets the initial uncertainty on the model fixed effects.
  void SetInitialUncertainty(const LinearVariableMapType &iu) { m_InitialUncertainty = iu; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Runs the gradient ascent algorithm and updates the statistical model.
  virtual void Update();

  /// Prints the algorithm current state.
  virtual void Print() const {
    Superclass::Print();
//    std::cout << "Log-likelihood = "
//              << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration].sum()
//              << "\t [data attachement = " << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration][0]
//              << " ; regularity = " << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration][1]
//              << "]" << std::endl;
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Shortcut function to compute the model functional. Sets the model fixed effects in the process.
  ScalarType ComputeFunctional(const VectorType &position);

  /// Golden ratio line search in the domain +/- direction. Returns the new functional.
  ScalarType LineSearchStep(const VectorType &direction,
                            const ScalarType functional,
                            VectorType &position);

  /// Expresses the initial uncertainty in convenient vector form.
  VectorType LinearizeInitialUncertainty(const LinearVariableMapType &modelFixedEffects) const;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Maximum number of iterations in the line search procedure.
  unsigned int m_MaxLineSearchIterations;

  /// Initial uncertainty on the fixed effects.
  LinearVariableMapType m_InitialUncertainty;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Working attribute(s), could be spared.
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Size parameters of the fixed effects, saved for unvectorization.
  std::vector<std::vector<unsigned int>> m_Structure;

  /// Keys of the fixed effects, saved for unvectorization.
  std::vector<std::string> m_Keys;

}; /* class PowellsMethod */


