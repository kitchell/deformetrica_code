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
#include "AbstractSampler.h"

/// Input-output files.
#include "MatrixDLM.h"

/// Support file
#include "LinearAlgebra.h"
#include "ProbabilityDistributions.h"

using namespace def::algebra;
using namespace def::proba;

/**
 *	\brief      McmcSaem object class
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    A model estimator is an algorithm which updates the fixed effects of a statistical model.
 */
template<class ScalarType, unsigned int Dimension>
class McmcSaem : public AbstractEstimator<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Model estimator type.
  typedef AbstractEstimator<ScalarType, Dimension> Superclass;

  /// Abstract statistical model type.
  typedef typename Superclass::StatisticalModelType StatisticalModelType;

  /// Deformable multi object type.
  typedef typename StatisticalModelType::DeformableMultiObjectType DeformableMultiObjectType;

  /// Abstract sampler type.
  typedef AbstractSampler<ScalarType, Dimension> AbstractSamplerType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Constructor.
  McmcSaem();

  /// Copy constructor.
  McmcSaem(const McmcSaem &other);

  /// Makes a copy of the object.
  virtual McmcSaem *Clone() { return new McmcSaem(*this); }

  /// Destructor.
  virtual ~McmcSaem();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Gets the sampler.
  std::shared_ptr<AbstractSamplerType> GetSampler() const { return m_Sampler; }
  /// Sets the sampler.
  void SetSampler(std::shared_ptr<AbstractSamplerType> s) { m_Sampler = s; }

  /// Sets the memory window size.
  void SetMemoryWindowSize(unsigned int const &n) { m_MemoryWindowSize = n; }

  /// Sets the adaptive flag to true.
  void SetAdaptiveMode() { m_AdaptiveMode = true; }
  /// Sets the adaptive flag to false.
  void UnsetAdaptiveMode() { m_AdaptiveMode = false; }
  /// Sets the memory window size.
  void SetAdaptMemoryWindowSize(unsigned int const &n) { m_AdaptMemoryWindowSize = n; }

  /// Sets the tempering flag to true.
  void SetUseTempering() { m_UseTempering = true; }
  /// Sets the tempering flag to false.
  void UnsetUseTempering() { m_UseTempering = false; }

  /// Sets the initial temprature.
  void SetInitialTemperature(const ScalarType &temperature) { m_InitialTemperature = temperature; }
  /// Sets the tempering duration ratio.
  void SetTemperingDurationRatio(const ScalarType &tdr) { m_TemperingDurationRatio = tdr; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Runs the MCMC SAEM algorithm and updates the statistical model.
  virtual void Update();

  /// Prints the algorithm current state.
  virtual void Print() const;

  /// Writes the model parameters trajectory.
  virtual void Write() const;

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initializes the number of burn-in iterations.
  void InitializeNumberOfBurnInIterations();

  /// Initializes the initial sufficient statistics.
  void InitializeSufficientStatistics();

  /// Initializes the acceptance rates memory and related attributes.
  void InitializeAcceptanceRatesInformation();
  /// Updates the acceptance rates memory.
  void UpdateAcceptanceRatesInformation();

  /// Computes the step size.
  ScalarType ComputeStepSize() const;

  /// Initializes the model parameters trajectory.
  void InitializeModelParametersTrajectory();
  /// Updates the model parameters trajectory.
  void UpdateModelParametersTrajectory();

  /// Initializes the temperature.
  void InitializeTemperature();
  /// Updates the temperature.
  void UpdateTemperature();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Protected attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sampler.
  std::shared_ptr<AbstractSamplerType> m_Sampler;

  /// Current sufficient statistics.
  LinearVariableMapType m_SufficientStatistics;

  /// Number of iterations without memory.
  unsigned int m_NbBurnInIterations;

  /// Names of the random effects. Population first, individual then, alphabetical ordering.
  std::vector<std::string> m_RandomEffectsNames;

  /// Acceptance rates of the current iteration.
  VectorType m_CurrentAcceptanceRates;
  /// Mean acceptance rates, computed over all past iterations.
  VectorType m_MeanAcceptanceRates;

  /// Size of the averaging window for the cumulated acceptance rates.
  unsigned int m_MemoryWindowSize;
  /// Memory of the last m_MemoryWindowSize acceptances rates.
  //  m_MemoryWindowSize rows x m_RandomEffectsNames.size() cols.
  MatrixType m_AcceptanceRatesMemory;
  /// Mean acceptance rates, computed over the last m_MemoryWindowSize iterations.
  VectorType m_WindowMeanAcceptanceRates;

  /// Flag that indicates if the sampler is to be dynamically adapted.
  bool m_AdaptiveMode;
  /// Size of the averaging window for the cumulated acceptance rates.
  unsigned int m_AdaptMemoryWindowSize;
  /// Memory of the last m_MemoryWindowSize acceptances rates.
  //  m_MemoryWindowSize rows x m_RandomEffectsNames.size() cols.
  MatrixType m_AcceptanceRatesAdaptMemory;
  /// Mean acceptance rates, computed over the last m_MemoryWindowSize iterations.
  VectorType m_AdaptWindowMeanAcceptanceRates;

  /// Memory of the model parameters along the estimation.
  MatrixType m_ModelParametersTrajectory;
  /// Resolution of the model parameters trajectory.
  unsigned int m_SaveModelParametersEveryNIters;

  /// Flag that indicates if tempering of model log-likelihood should be used.
  bool m_UseTempering;
  /// Current temperature parameter.
  ScalarType m_CurrentTemperature;
  /// Initial temperature parameter.
  ScalarType m_InitialTemperature;
  /// Number of iterations with a temperature greater than one, expressed as a ratio of m_NbBurnInIterations.
  ScalarType m_TemperingDurationRatio;
  /// /// Number of iterations with a temperature greater than one.
  unsigned int m_NbTemperingIterations;
  /// Temperature update parameter.
  ScalarType m_TemperatureUpdateParameter;

}; /* class McmcSaem */


