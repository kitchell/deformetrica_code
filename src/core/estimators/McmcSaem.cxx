/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "McmcSaem.h"
#include "MatrixDLM.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
McmcSaem<ScalarType, Dimension>
::McmcSaem() : Superclass(), m_NbBurnInIterations(100), m_MemoryWindowSize(100),
               m_UseTempering(false), m_InitialTemperature(100.0), m_TemperingDurationRatio(0.5) {
  this->SetMcmcSaemType();
}

template<class ScalarType, unsigned int Dimension>
McmcSaem<ScalarType, Dimension>
::~McmcSaem() {}

template<class ScalarType, unsigned int Dimension>
McmcSaem<ScalarType, Dimension>
::McmcSaem(const McmcSaem &other) {
  m_Sampler = other.m_Sampler;
  m_SufficientStatistics = other.m_SufficientStatistics;
  m_NbBurnInIterations = other.m_NbBurnInIterations;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::Update() {

  /// Declare variables.
  unsigned int iter(1);
  ScalarType step;
  LinearVariableMapType averagedPopRER, sufficientStatistics, fixedEffects;
  LinearVariablesMapType averagedIndRER;

  /// Auxiliary function.
  auto check_file = [](const std::string &f) {
    std::ifstream infile(f);
    return infile.good();
  };

  /// Initialization.
  InitializeNumberOfBurnInIterations();
  InitializeAcceptanceRatesInformation();
  InitializeSufficientStatistics();
  InitializeModelParametersTrajectory();
  InitializeTemperature();

  bool computation_end_state = false;
  def::utils::DeformationState state;
  if (def::utils::settings.load_state && check_file(def::utils::settings.input_state_filename)) {

    /* SERIALIZE : LOAD */
    state.load(def::utils::settings.input_state_filename);
    state >> computation_end_state;

    if (computation_end_state) {
      std::cout << std::endl << "WARNING: Computation was already completed using the current serialize file"
                << std::endl;
    }

    state >> iter >> Superclass::m_PopulationRER >> Superclass::m_IndividualRER
          >> averagedPopRER >> averagedIndRER
          >> m_SufficientStatistics
          >> m_CurrentAcceptanceRates >> m_MeanAcceptanceRates
          >> m_AcceptanceRatesMemory >> m_WindowMeanAcceptanceRates
          >> m_AcceptanceRatesAdaptMemory >> m_AdaptWindowMeanAcceptanceRates
          >> m_ModelParametersTrajectory
          >> m_CurrentTemperature >> m_NbTemperingIterations >> m_TemperatureUpdateParameter
          >> fixedEffects;

    Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);
    m_Sampler->RecoverState(state);

    /// Print console information.
    Superclass::m_CurrentIteration = iter;
    std::cout << "\n[ Loading saved state : iteration " << iter << " ]" << std::endl;
    Print();

    /// Sets the iteration number.
    iter += 1;
    Superclass::m_CurrentIteration = iter;

  } else {

    /* STANDARD INITIALIZATION */
    /// Ensures that all the model fixed effects are initialized.
    VectorType logLikelihoodTerms;
    Superclass::m_StatisticalModel
        ->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
            Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER, logLikelihoodTerms);

    /// Print console information.
    Superclass::Print();
    std::cout << "Log-likelihood = " << logLikelihoodTerms.sum() << "\t (+ cst)" << std::endl;

    std::cout << ">> MCMC-SAEM algorithm launched for " << Superclass::m_MaxIterations
              << " iterations (" << m_NbBurnInIterations << " iterations of burn-in)" << std::endl;

    if (m_UseTempering) {
      std::cout << ">> Tempering option activated :" << std::endl;
      std::cout << "\t\tThe optimized log-likelihood will be tempered during the first "
                << m_NbTemperingIterations << " iterations" << std::endl;
      std::cout << "\t\tInitialization of the temperature : T = " << m_CurrentTemperature << std::endl;
      std::cout << "\t\tInitialization of the temperature update parameter : zeta = "
                << m_TemperatureUpdateParameter << std::endl;
    }

    /// Initialization of the average random effects realizations.
    averagedPopRER = Superclass::m_PopulationRER;
    averagedPopRER.fill(0.0);
    averagedIndRER = Superclass::m_IndividualRER;
    averagedPopRER.fill(0.0);
  }

  m_Sampler->Initialize(Superclass::m_PopulationRER, Superclass::m_IndividualRER);

  /// Main loop.
  for (; iter < Superclass::m_MaxIterations + 1; ++iter) {
    Superclass::m_CurrentIteration = iter;

    /// Simulation.
    m_Sampler->Sample(Superclass::m_PopulationRER, Superclass::m_IndividualRER,
                      m_CurrentAcceptanceRates, m_CurrentTemperature);

    /// Stochastic approximation.
    Superclass::m_StatisticalModel->ComputeSufficientStatistics(Superclass::m_DataSet,
                                                                Superclass::m_PopulationRER,
                                                                Superclass::m_IndividualRER,
                                                                sufficientStatistics);

    step = ComputeStepSize();
    m_SufficientStatistics += step * (sufficientStatistics - m_SufficientStatistics);

    /// Maximization and update.
    Superclass::m_StatisticalModel->UpdateFixedEffects(
        Superclass::m_DataSet, m_SufficientStatistics, m_CurrentTemperature);

    /// Averages the random effect realizations in the concentration phase.
    if (step < 1.0) {
      averagedPopRER *= (iter - m_NbBurnInIterations) / (iter + 1 - m_NbBurnInIterations);
      averagedPopRER += Superclass::m_PopulationRER / (iter + 1 - m_NbBurnInIterations);
      averagedIndRER *= (iter - m_NbBurnInIterations) / (iter + 1 - m_NbBurnInIterations);
      averagedIndRER += Superclass::m_IndividualRER / (iter + 1 - m_NbBurnInIterations);
    }

    /// Displays information about the current state of the algorithm.
    UpdateAcceptanceRatesInformation();
    if (!(iter % m_SaveModelParametersEveryNIters)) { UpdateModelParametersTrajectory(); }
    if (!(iter % Superclass::m_PrintEveryNIters)) { Print(); } // Prints information.
    if (!(iter % Superclass::m_SaveEveryNIters)) {             // Saves information.
      /* SERIALIZATION */
      if (def::utils::settings.save_state) {
        Superclass::m_StatisticalModel->GetFixedEffects(fixedEffects);
        state << computation_end_state << iter << Superclass::m_PopulationRER << Superclass::m_IndividualRER
              << averagedPopRER << averagedIndRER
              << m_SufficientStatistics
              << m_CurrentAcceptanceRates << m_MeanAcceptanceRates
              << m_AcceptanceRatesMemory << m_WindowMeanAcceptanceRates
              << m_AcceptanceRatesAdaptMemory << m_AdaptWindowMeanAcceptanceRates
              << m_ModelParametersTrajectory
              << m_CurrentTemperature << m_NbTemperingIterations << m_TemperatureUpdateParameter
              << fixedEffects;
        m_Sampler->SaveState(state);
        state.save_and_reset(def::utils::settings.output_state_filename);
      }

      Superclass::m_StatisticalModel->Write(Superclass::m_DataSet,
                                            Superclass::m_PopulationRER, Superclass::m_IndividualRER);
      Write();
    }
    if (m_AdaptiveMode && !(iter % m_AdaptMemoryWindowSize)) { // Optionally adapts the sampler.
      m_Sampler->AdaptProposalDistributions(
          m_AdaptWindowMeanAcceptanceRates, iter, !(iter % Superclass::m_PrintEveryNIters));
    }
    if (m_UseTempering) { UpdateTemperature(); }
  }

  /// Write outputs.
  std::cout << "Write output files ...";

  /* SERIALIZATION */
  computation_end_state = true;
  if (def::utils::settings.save_state) {
    Superclass::m_StatisticalModel->GetFixedEffects(fixedEffects);
    state << computation_end_state << iter - 1 << Superclass::m_PopulationRER << Superclass::m_IndividualRER
          << averagedPopRER << averagedIndRER
          << m_SufficientStatistics
          << m_CurrentAcceptanceRates << m_MeanAcceptanceRates
          << m_AcceptanceRatesMemory << m_WindowMeanAcceptanceRates
          << m_AcceptanceRatesAdaptMemory << m_AdaptWindowMeanAcceptanceRates
          << m_ModelParametersTrajectory
          << m_CurrentTemperature << m_NbTemperingIterations << m_TemperatureUpdateParameter
          << fixedEffects;
    m_Sampler->SaveState(state);
    state.save_and_reset(def::utils::settings.output_state_filename);
  }

  Superclass::m_StatisticalModel->Write(Superclass::m_DataSet, averagedPopRER, averagedIndRER);
  Write();

  std::cout << "... done." << std::endl;
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::Print() const {
  /// Prints the iteration number.
  Superclass::Print();

  /// if the tempering option is active, prints the current temperature.
  std::cout << ">> Current temperature : T = " << m_CurrentTemperature << std::endl;

  /// Initialization.
  unsigned int k = 0;

//  /// Prints the instantaneous acceptance rates.
//  if (Superclass::m_PrintEveryNIters == 1) {
//    std::cout << ">> Instantaneous acceptance rates = ";
//    for (auto it = m_RandomEffectsNames.begin(); it != m_RandomEffectsNames.end(); ++it, ++k) {
//      std::cout << std::endl;
//      std::cout << "\t\t" << m_CurrentAcceptanceRates(k) << " % \t[" << *it << "]";
//    }
//    std::cout << std::endl;
//  }
//
//  if (m_AdaptiveMode && m_AdaptMemoryWindowSize != 1) {
//    /// Prints the mean acceptance rates, averaged over the m_AdaptMemoryWindowSize last iterations.
//    k = 0;
//    std::cout << ">> Mean acceptance rates (last " << m_AdaptMemoryWindowSize << " iterations) = ";
//    for (auto it = m_RandomEffectsNames.begin(); it != m_RandomEffectsNames.end(); ++it, ++k) {
//      std::cout << std::endl;
//      std::cout << "\t\t" << m_AdaptWindowMeanAcceptanceRates(k) << " % \t[" << *it << "]";
//    }
//    std::cout << std::endl;
//  }
//
//  /// Prints the mean acceptance rates, averaged over the m_MemoryWindowSize last iterations.
//  k = 0;
//  std::cout << ">> Mean acceptance rates (last " << m_MemoryWindowSize << " iterations) = ";
//  for (auto it = m_RandomEffectsNames.begin(); it != m_RandomEffectsNames.end(); ++it, ++k) {
//    std::cout << std::endl;
//    std::cout << "\t\t" << m_WindowMeanAcceptanceRates(k) << " % \t[" << *it << "]";
//  }
//  std::cout << std::endl;

  /// Prints the mean acceptance rates, averaged over all the past iterations.
  k = 0;
  std::cout << ">> Mean acceptance rates (all past iterations) = ";
  for (auto it = m_RandomEffectsNames.begin(); it != m_RandomEffectsNames.end(); ++it, ++k) {
    std::cout << std::endl;
    std::cout << "\t\t" << m_MeanAcceptanceRates(k) << " % \t[" << *it << "]";
  }
  std::cout << std::endl;

  /// Prints information about the current statistical model.
  Superclass::m_StatisticalModel->Print();
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::Write() const {
  const std::string outputDir = def::utils::settings.output_dir;

  std::ostringstream oss;
  oss << outputDir << Superclass::m_StatisticalModel->GetName() << "__ParametersTrajectory.txt" << std::ends;

  writeMatrixDLM<ScalarType>(oss.str().c_str(), m_ModelParametersTrajectory.get_n_rows(
      0, Superclass::m_CurrentIteration / m_SaveModelParametersEveryNIters + 1));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::InitializeNumberOfBurnInIterations() {
  if (Superclass::m_MaxIterations > 4000) { m_NbBurnInIterations = Superclass::m_MaxIterations - 2000; }
  else { m_NbBurnInIterations = Superclass::m_MaxIterations / 2; }
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::InitializeSufficientStatistics() {
  Superclass::m_StatisticalModel->ComputeSufficientStatistics(Superclass::m_DataSet,
                                                              Superclass::m_PopulationRER,
                                                              Superclass::m_IndividualRER,
                                                              m_SufficientStatistics);
  m_SufficientStatistics.fill(0.0);
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::InitializeAcceptanceRatesInformation() {
  /// Names of the random effects.
  unsigned int k = 0;
  for (auto it = Superclass::m_PopulationRER.begin(); it != Superclass::m_PopulationRER.end(); ++it, ++k) {
    m_RandomEffectsNames.push_back(it->first);
  }
  for (auto it = Superclass::m_IndividualRER.begin(); it != Superclass::m_IndividualRER.end(); ++it, ++k) {
    m_RandomEffectsNames.push_back(it->first);
  }

  /// Initializes the current acceptance rates vector.
  m_CurrentAcceptanceRates.set_size(k);
  m_CurrentAcceptanceRates.fill(0.0);

  /// Initializes the acceptance rates memory matrix.
  m_AcceptanceRatesMemory.set_size(m_MemoryWindowSize, k);
  m_AcceptanceRatesMemory.fill(0.0);
  /// Initializes the window mean acceptance rates vector.
  m_WindowMeanAcceptanceRates.set_size(k);
  m_WindowMeanAcceptanceRates.fill(0.0);

  if (m_AdaptiveMode) {
    /// Initializes the acceptance rates adapt memory matrix.
    m_AcceptanceRatesAdaptMemory.set_size(m_AdaptMemoryWindowSize, k);
    m_AcceptanceRatesAdaptMemory.fill(0.0);
    /// Initializes the adapt window mean acceptance rates vector.
    m_AdaptWindowMeanAcceptanceRates.set_size(k);
    m_AdaptWindowMeanAcceptanceRates.fill(0.0);
  }

  /// Initializes the mean acceptance rates vector.
  m_MeanAcceptanceRates.set_size(k);
  m_MeanAcceptanceRates.fill(0.0);
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::UpdateAcceptanceRatesInformation() {
  /// Updates the memory of the last m_MemoryWindowSize acceptance rates.
  m_AcceptanceRatesMemory.set_row((Superclass::m_CurrentIteration - 1) % m_MemoryWindowSize,
                                  m_CurrentAcceptanceRates);

  /// Updates the window mean acceptance rates vector.
  for (unsigned int k = 0; k < m_RandomEffectsNames.size(); ++k)
    m_WindowMeanAcceptanceRates(k) = m_AcceptanceRatesMemory.get_column(k).sum()
        / std::min(m_MemoryWindowSize, Superclass::m_CurrentIteration);

  if (m_AdaptiveMode) {
    /// Updates the adapt memory of the last m_AdaptMemoryWindowSize acceptance rates.
    m_AcceptanceRatesAdaptMemory.set_row((Superclass::m_CurrentIteration - 1) % m_AdaptMemoryWindowSize,
                                         m_CurrentAcceptanceRates);

    /// Updates the window mean acceptance rates vector.
    for (unsigned int k = 0; k < m_RandomEffectsNames.size(); ++k)
      m_AdaptWindowMeanAcceptanceRates(k) = m_AcceptanceRatesAdaptMemory.get_column(k).sum()
          / std::min(m_AdaptMemoryWindowSize, Superclass::m_CurrentIteration);
  }

  /// Updates the mean acceptance rates vector.
  m_MeanAcceptanceRates *= 1.0 - 1.0 / Superclass::m_CurrentIteration;
  m_MeanAcceptanceRates += m_CurrentAcceptanceRates / Superclass::m_CurrentIteration;
}

template<class ScalarType, unsigned int Dimension>
ScalarType
McmcSaem<ScalarType, Dimension>
::ComputeStepSize() const {
  int aux = Superclass::m_CurrentIteration - m_NbBurnInIterations + 1;
  if (aux <= 0) { return 1.0; }
  else { return pow(aux, -0.6); }
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::InitializeModelParametersTrajectory() {
  const unsigned int totalNumberOfTrajectoryPoints = 1000;
  m_SaveModelParametersEveryNIters = Superclass::m_MaxIterations / totalNumberOfTrajectoryPoints;
  if (m_SaveModelParametersEveryNIters == 0) { m_SaveModelParametersEveryNIters = 1; }

  LinearVariableMapType modelParameters;
  Superclass::m_StatisticalModel->GetFixedEffects(modelParameters);

  VectorType vectorizedModelParameters(0);
  for (auto it = modelParameters.begin(); it != modelParameters.end(); ++it) {
    vectorizedModelParameters.push_back(it->second.vectorize());
  }

  m_ModelParametersTrajectory.set_size(totalNumberOfTrajectoryPoints + 1, vectorizedModelParameters.size());
  m_ModelParametersTrajectory.set_row(0, vectorizedModelParameters);
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::UpdateModelParametersTrajectory() {
  LinearVariableMapType modelParameters;
  Superclass::m_StatisticalModel->GetFixedEffects(modelParameters);

  VectorType vectorizedModelParameters(0);
  for (auto it = modelParameters.begin(); it != modelParameters.end(); ++it) {
    vectorizedModelParameters.push_back(it->second.vectorize());
  }

  m_ModelParametersTrajectory.set_row(Superclass::m_CurrentIteration / m_SaveModelParametersEveryNIters,
                                      vectorizedModelParameters);
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::InitializeTemperature() {
  if (m_UseTempering) {
    assert(m_InitialTemperature > 1);
    m_CurrentTemperature = m_InitialTemperature;

    assert (m_TemperingDurationRatio > 0 && m_TemperingDurationRatio < 1);
    m_NbTemperingIterations = m_TemperingDurationRatio * (m_NbBurnInIterations + 1);
    m_TemperatureUpdateParameter = pow(m_InitialTemperature, -1.0 / (1e-5 + m_NbTemperingIterations));
  } else {
    m_CurrentTemperature = 1.0;
  }
}

template<class ScalarType, unsigned int Dimension>
void
McmcSaem<ScalarType, Dimension>
::UpdateTemperature() {
  if (Superclass::m_CurrentIteration <= m_NbTemperingIterations) {
    m_CurrentTemperature *= m_TemperatureUpdateParameter;
    assert(m_CurrentTemperature > 1.0);
  } else { m_CurrentTemperature = 1.0; }
}

template
class McmcSaem<double, 2>;
template
class McmcSaem<double, 3>;

