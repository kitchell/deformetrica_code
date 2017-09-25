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
#include "ParametricImage.h"
#include "LinearAlgebra.h"

using namespace def::algebra;

/**
 *	\brief      GradientAscent object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A model estimator is an algorithm which updates the fixed effects of a statistical model.
 */
template<class ScalarType, unsigned int Dimension>
class GradientAscent : public AbstractEstimator<ScalarType, Dimension> {
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
  GradientAscent();

  /// Copy constructor.
  GradientAscent(const GradientAscent &other);

  /// Makes a copy of the object.
  virtual GradientAscent *Clone() { return new GradientAscent(*this); }

  /// Destructor.
  virtual ~GradientAscent();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns the maximum number of iterations in the line search procedure.
  inline unsigned int GetMaxLineSearchIterations() const { return m_MaxLineSearchIterations; }
  /// Sets the maximum number of iterations in the line search procedure to \e n.
  inline void SetMaxLineSearchIterations(const unsigned int n) { m_MaxLineSearchIterations = n; }

  /// Returns the initial step size.
  inline ScalarType GetInitialStepSize() const { return m_InitialStepSize; }
  /// Sets the initial step size to \e x.
  inline void SetInitialStepSize(const ScalarType x) { m_InitialStepSize = x; }

  /// Returns the shrink parameter for the line search.
  inline ScalarType GetAdaptiveShrink() const { return m_AdaptiveShrink; }
  /// Sets the shrink parameter for the line search \e x.
  inline void SetAdaptiveShrink(const ScalarType x) { m_AdaptiveShrink = x; }

  /// Returns the expand parameter for the line search.
  inline ScalarType GetAdaptiveExpand() const { return m_AdaptiveExpand; }
  /// Sets the expand parameter for the line search to \e x.
  inline void SetAdaptiveExpand(const ScalarType x) { m_AdaptiveExpand = x; }

  /// Returns the adaptive tolerance parameter.
  inline ScalarType GetAdaptiveTolerance() const { return m_AdaptiveTolerance; }
  /// Sets the adaptive tolerance parameter to \e x.
  inline void SetAdaptiveTolerance(const ScalarType x) { m_AdaptiveTolerance = x; }

  /// Sets the multivariate line search flag.
  void SetMultivariateLineSearchFlag(bool const &flag) { m_MultivariateLineSearchFlag = flag; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Runs the gradient ascent algorithm and updates the statistical model.
  virtual void Update();

  /// Prints the algorithm current state.
  virtual void Print() const {
    Superclass::Print();
    std::cout << "Log-likelihood = "
              << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration].sum()
              << "\t [data attachement = " << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration][0]
              << " ; regularity = " << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration][1]
              << "]" << std::endl;
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
  * \brief		Computes a gradient ascent step for test values of the variables and step sizes
  *
  * \param[out]	newModelParameters	  Updated value of the model parameters.
  * \param[in]	modelParameters		  Current value of the model parameters.
  * \param[out]	functionalGradient	  Gradient of the functional wrt the parameters.
  * \param[in]	step			      Step-sizes for parameters update.
  */
  void GradientAscentStep(const LinearVariableMapType &fixedEffects,
                          const LinearVariableMapType &popGrad,
                          const LinearVariablesMapType &indGrad,
                          const VectorType &step,
                          LinearVariableMapType &newFixedEffects,
                          LinearVariableMapType &newPopRER,
                          LinearVariablesMapType &newIndRER);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Vector containing a history of the log-likelihood values during the optimization method.
  std::vector<VectorType> m_LogLikelihoodTermsHistory;

  /// Maximum number of iterations in the line search procedure.
  unsigned int m_MaxLineSearchIterations;

  /// Initial step size.
  ScalarType m_InitialStepSize;

  /// Shrink parameter for the line search.
  ScalarType m_AdaptiveShrink;
  /// Expand parameter for the line search.
  ScalarType m_AdaptiveExpand;
  /// The algorithm stops when \f$F(end-1)-F(end) < \verb#m_AdaptiveTolerance# * \left( F(0)-F(end) \right)\f$.
  ScalarType m_AdaptiveTolerance;

  /// Flag that triggers a multivariate line search.
  bool m_MultivariateLineSearchFlag;

}; /* class GradientAscent */


