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

/// Core files.
#include "AbstractStatisticalModel.h"
#include "LongitudinalDataSet.h"

/// Support files.
#include "LinearAlgebra.h"
#include <src/support/utilities/GeneralSettings.h>
#include <src/support/utilities/SerializeDeformationState.h>

/// Standard files.
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

using namespace def::algebra;

/**
 *	\brief      AbstractModelEstimator object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A model estimator is an algorithm which updates the fixed effects of a statistical model.
 */
template<class ScalarType, unsigned int Dimension>
class AbstractEstimator {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Typedefs :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Possible types of model estimators.
  typedef enum {
    null,               // Null value.
    GradientAscent,     // Gradient ascent.
    FastGradientAscent, // Fast gradient ascent.
    McmcSaem,           // Stochastic approximation expectancy maximization.
    PowellsMethod        // Powell's method.
  } ModelEstimatorType;

  /// Abstract statistical model type.
  typedef AbstractStatisticalModel<ScalarType, Dimension> StatisticalModelType;

  /// Longitudinal dataset type.
  typedef LongitudinalDataSet<ScalarType, Dimension> LongitudinalDataSetType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  AbstractEstimator();

  /// Copy constructor.
  AbstractEstimator(const AbstractEstimator &other);

  /// Makes a copy of the object.
  virtual AbstractEstimator *Clone() = 0;

  /// Destructor
  virtual ~AbstractEstimator();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Returns true if the model estimator is a gradient ascent.
  bool IsGradientAscent() const { return (m_Type == GradientAscent); }
  /// Sets gradient ascent type.
  void SetGradientAscentType() { m_Type = GradientAscent; }

  /// Returns true if the model estimator is a fast gradient ascent.
  bool IsFastGradientAscent() const { return (m_Type == FastGradientAscent); }
  /// Sets gradient ascent type.
  void SetFastGradientAscentType() { m_Type = FastGradientAscent; }

  /// Returns true if the model estimator is a fast gradient ascent.
  bool IsMcmcSaem() const { return (m_Type == McmcSaem); }
  /// Sets gradient ascent type.
  void SetMcmcSaemType() { m_Type = McmcSaem; }

  /// Returns true if the model estimator is a Powell's algorithm.
  bool IsPowellsMethod() const { return (m_Type == PowellsMethod); }
  /// Sets Powell's method type.
  void SetPowellsMethodType() { m_Type = PowellsMethod; }

  /// Returns the statistical model.
  std::shared_ptr<StatisticalModelType> GetStatisticalModel() const { return m_StatisticalModel; }
  /// Sets the statistical model to \e model.
  void SetStatisticalModel(std::shared_ptr<StatisticalModelType> const model) { m_StatisticalModel = model; }

  /// Returns the data set.
  //StatisticalModelType* GetDataSet() const { return m_DataSet; }
  /// Sets the data set to \e dataSet
  void SetDataSet(LongitudinalDataSetType *dataSet) {
    m_DataSet = dataSet;
    m_NumberOfSubjects = dataSet->GetNumberOfSubjects();
  }

  /// Returns the maximum of iterations.
  unsigned int GetMaxIterations() const { return m_MaxIterations; }
  /// Sets the maximum of iterations to \e n.
  void SetMaxIterations(const unsigned int n) { m_MaxIterations = n; }

  /// Returns the number of iterations between two savings.
  unsigned int GetSaveEveryNIters() const { return m_SaveEveryNIters; }
  /// Sets the number of iterations between two savings to \e n
  void SetSaveEveryNIters(const unsigned int n) { m_SaveEveryNIters = n; }

  /// Sets the number of iterations between two printing to \e n
  void SetPrintEveryNIters(const unsigned int n) { m_PrintEveryNIters = n; }

  /// Sets the initial population random effects realization ("RER").
  void InitializePopulationRER(LinearVariableMapType const &map) { m_PopulationRER = map; }
  /// Sets the initial individual random effects realization ("RER").
  void InitializeIndividualRER(LinearVariablesMapType const &map) { m_IndividualRER = map; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Runs the algorithm and updates the statistical model.
  virtual void Update() = 0;

  /// Prints the algorithm current state.
  virtual void Print() const {
    std::cout << "\n------------------------------------- Iteration " << m_CurrentIteration
              << " -------------------------------------" << std::endl;
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Type of the model estimator.
  ModelEstimatorType m_Type;

  /// Target statistical model.
  std::shared_ptr<StatisticalModelType> m_StatisticalModel;

  /// Target data set.
  LongitudinalDataSetType *m_DataSet;
  /// Number of subjects in the dataset.
  unsigned int m_NumberOfSubjects;

  /// Maximum number of iterations.
  unsigned int m_MaxIterations;
  /// The current model state will be saved every m_SaveEveryNIters iterations.
  unsigned int m_SaveEveryNIters;
  /// The algorithm will print relevant information every m_PrintEveryNIters iterations.
  unsigned int m_PrintEveryNIters;

  /// Current iteration.
  unsigned int m_CurrentIteration;

  /// Current population random effects realization ("RER").
  LinearVariableMapType m_PopulationRER;
  /// Current individual random effects realization ("RER").
  LinearVariablesMapType m_IndividualRER;

}; /* class AbstractEstimator */

