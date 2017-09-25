/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/



#include "FastGradientAscent.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
FastGradientAscent<ScalarType, Dimension>
::FastGradientAscent() : Superclass(), m_NbIterFreeze(3), m_InitialStepSize(0.001), m_MultivariateLineSearchFlag(true) {
  this->SetFastGradientAscentType();
}

template<class ScalarType, unsigned int Dimension>
FastGradientAscent<ScalarType, Dimension>
::~FastGradientAscent() {}

template<class ScalarType, unsigned int Dimension>
FastGradientAscent<ScalarType, Dimension>
::FastGradientAscent(const FastGradientAscent &other) {
  m_MaxLineSearchIterations = other.m_MaxLineSearchIterations;
  m_AdaptiveShrink = other.m_AdaptiveShrink;
  m_AdaptiveExpand = other.m_AdaptiveExpand;
  m_AdaptiveTolerance = other.m_AdaptiveTolerance;
  m_InitialStepSize = other.m_InitialStepSize;
  m_MultivariateLineSearchFlag = other.m_MultivariateLineSearchFlag;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
FastGradientAscent<ScalarType, Dimension>
::Update() {

  /// Declare variables.
  LinearVariableMapType fixedEffects, auxFixedEffects, newFixedEffects;
  LinearVariableMapType newPopRER(Superclass::m_PopulationRER), auxPopRER(Superclass::m_PopulationRER);
  LinearVariablesMapType newIndRER(Superclass::m_IndividualRER), auxIndRER(Superclass::m_IndividualRER);

  LinearVariableMapType popGrad;
  LinearVariablesMapType indGrad;
  VectorType gradSqNorms;

  unsigned int iter(1), iterRef(0), k(0), freezeDirectionCounter(1);
  unsigned long nbParams;
  ScalarType lsqRef, tau(1.0);
  VectorType newLogLikelihoodTerms, step;

  /// Auxiliary function.
  auto check_file = [](const std::string &f) {
    std::ifstream infile(f);
    return infile.good();
  };

  /// Initialization.
  bool computation_end_state = false;
  def::utils::DeformationState deformation_state;
  if (def::utils::settings.load_state && check_file(def::utils::settings.input_state_filename)) {

    /* SERIALIZE : LOAD */
    deformation_state.load(def::utils::settings.input_state_filename);
    deformation_state >> computation_end_state;

    if (computation_end_state)
      std::cout << std::endl << "WARNING: Computation was already completed using the current serialize file" << std::endl;
    else
      std::cout << std::endl << "Loading serialization file..." << std::endl;


    deformation_state >> iter >> freezeDirectionCounter >> tau
                      >> step >> gradSqNorms >> m_LogLikelihoodTermsHistory
                      >> fixedEffects >> auxFixedEffects
                      >> Superclass::m_PopulationRER >> auxPopRER
                      >> Superclass::m_IndividualRER >> auxIndRER
                      >> popGrad >> indGrad;

    Superclass::m_StatisticalModel->SetFixedEffects(auxFixedEffects);
    nbParams = popGrad.size() + indGrad.size();
    lsqRef = m_LogLikelihoodTermsHistory[iter].sum();

    std::cout << "\n--------------------------------- Loading iteration " << iter
              << " ---------------------------------" << std::endl;
    std::cout << "Log-likelihood = "
              << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration].sum()
              << "\t [data attachement = " << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration][0]
              << " ; regularity = " << m_LogLikelihoodTermsHistory[Superclass::m_CurrentIteration][1]
              << "]" << std::endl;

    iter += 1;
    Superclass::m_CurrentIteration = iter;

  } else {

    /* STANDARD INITIALIZATION */
    Superclass::m_StatisticalModel->GetFixedEffects(fixedEffects);
    auxFixedEffects = fixedEffects;

    m_LogLikelihoodTermsHistory.resize(Superclass::m_MaxIterations + 1);
    Superclass::m_StatisticalModel->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
        Superclass::m_DataSet, auxPopRER, auxIndRER, m_LogLikelihoodTermsHistory[0]);
    lsqRef = m_LogLikelihoodTermsHistory[0].sum();

    Print();

    Superclass::m_StatisticalModel
        ->ComputeCompleteLogLikelihoodGradient(Superclass::m_DataSet, auxPopRER, auxIndRER,
                                               popGrad, indGrad, gradSqNorms);
    nbParams = popGrad.size() + indGrad.size();
    step.set_size(nbParams);
    step.fill(m_InitialStepSize);
  }

  /// Main loop.
  for (; iter < Superclass::m_MaxIterations + 1; iter++) {
    Superclass::m_CurrentIteration = iter;

    bool foundMin = false;
    for (unsigned int li = 0; li < m_MaxLineSearchIterations; li++) {
      if (!(iter % Superclass::m_PrintEveryNIters)) {
        k = 0;
        std::cout << "Step size  = ";
        for (auto it = popGrad.begin(); it != popGrad.end(); ++it, ++k)
          std::cout << "\t " << step(k) << " [" << it->first << "] ";
        for (auto it = indGrad.begin(); it != indGrad.end(); ++it, ++k)
          std::cout << "\t " << step(k) << " [" << it->first << "] ";
        std::cout << std::endl;
      }

      GradientAscentStep(auxFixedEffects, auxPopRER, auxIndRER, popGrad, indGrad, step,
                         newFixedEffects, newPopRER, newIndRER);

      Superclass::m_StatisticalModel->SetFixedEffects(newFixedEffects);
      bool oob = Superclass::m_StatisticalModel
          ->UpdateFixedEffectsAndComputeCompleteLogLikelihood(Superclass::m_DataSet,
                                                              newPopRER,
                                                              newIndRER,
                                                              newLogLikelihoodTerms);

      ScalarType Q = newLogLikelihoodTerms.sum() - lsqRef - dot_product(gradSqNorms, step) * m_AdaptiveShrink;

      if (Q >= 0 && !oob) // Easy case : J < Jcurr.
      {
        foundMin = true;
        if (li == 0)
          freezeDirectionCounter++;
        break;
      } else if (nbParams > 1 && m_MultivariateLineSearchFlag) {
        // All but one step are decreased. All configurations are tried.
        freezeDirectionCounter = 1;
        step *= m_AdaptiveShrink;
        std::vector<LinearVariableMapType> newFixedEffects_prop(nbParams);
        std::vector<LinearVariableMapType> newPopRER_prop(nbParams);
        std::vector<LinearVariablesMapType> newIndRER_prop(nbParams);
        std::vector<VectorType> logLikelihoodTerms_prop(nbParams);
        std::vector<bool> oob_prop(nbParams);
        VectorType Q(nbParams);

        for (unsigned int k = 0; k < nbParams; ++k) {
          VectorType localStep(step);
          localStep(k) /= m_AdaptiveShrink;

          GradientAscentStep(auxFixedEffects, auxPopRER, auxIndRER, popGrad, indGrad, localStep,
                             newFixedEffects_prop[k], newPopRER_prop[k], newIndRER_prop[k]);

          Superclass::m_StatisticalModel->SetFixedEffects(newFixedEffects_prop[k]);
          oob_prop[k] = Superclass::m_StatisticalModel->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
              Superclass::m_DataSet, newPopRER_prop[k], newIndRER_prop[k], logLikelihoodTerms_prop[k]);

          Q[k] = logLikelihoodTerms_prop[k].sum() - lsqRef - dot_product(gradSqNorms, localStep) * m_AdaptiveShrink;
        }

        unsigned int index = Q.index_max();

        if (Q[index] >= 0 && !oob_prop[index]) {
          newLogLikelihoodTerms = logLikelihoodTerms_prop[index];
          newFixedEffects = newFixedEffects_prop[index];
          newPopRER = newPopRER_prop[index];
          newIndRER = newIndRER_prop[index];
          step(index) /= m_AdaptiveShrink;
          foundMin = true;
          break;
        }
      } else {
        freezeDirectionCounter = 1;
        step *= m_AdaptiveShrink;
      }
    } // Line search loop.

    if (!foundMin) // Line search loop terminated without finding a larger log-likelihood.
    {
      Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);
      std::cout << "Number of line search loops exceeded." << std::endl;
      break;
    }

    m_LogLikelihoodTermsHistory[iter] = newLogLikelihoodTerms;

    ScalarType tau_next = (1 + sqrt(1.0 + 4.0 * tau * tau)) / 2.0;
    if (freezeDirectionCounter > m_NbIterFreeze) // Check if we're in an increasing trend.
    {
      step *= m_AdaptiveExpand;
      tau_next = (1 + sqrt(1.0 + 4.0 * tau * tau / m_AdaptiveExpand)) / 2.0;
    }

    ScalarType tau_scale = (tau - 1.0) / tau_next;

    for (auto it = popGrad.begin(); it != popGrad.end(); ++it, ++k) {
      if (auxFixedEffects.count(it->first)) {
        auxFixedEffects[it->first]
            = newFixedEffects[it->first]
            + (newFixedEffects[it->first] - fixedEffects[it->first]) * tau_scale;
      } else {
        auxPopRER[it->first]
            = newPopRER[it->first]
            + (newPopRER[it->first] - Superclass::m_PopulationRER[it->first]) * tau_scale;
      }
    }
    for (auto it = indGrad.begin(); it != indGrad.end(); ++it, ++k) {
      for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i) {
        auxIndRER[it->first][i]
            = newIndRER[it->first][i]
            + (newIndRER[it->first][i] - Superclass::m_IndividualRER[it->first][i]) * tau_scale;
      }
    }

    fixedEffects = newFixedEffects;
    Superclass::m_PopulationRER = newPopRER;
    Superclass::m_IndividualRER = newIndRER;

    tau = tau_next;

    // Update the log-likelihood.
    Superclass::m_StatisticalModel->SetFixedEffects(auxFixedEffects);
    bool oob = Superclass::m_StatisticalModel->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
        Superclass::m_DataSet, auxPopRER, auxIndRER, newLogLikelihoodTerms);
    lsqRef = newLogLikelihoodTerms.sum();

    /// Displays information about the current state of the algorithm.
    if (!(iter % Superclass::m_PrintEveryNIters)) { Print(); }
    if (!(iter % Superclass::m_SaveEveryNIters)) {
      Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);
      Superclass::m_StatisticalModel->Write(Superclass::m_DataSet,
                                            Superclass::m_PopulationRER, Superclass::m_IndividualRER);
      Superclass::m_StatisticalModel->SetFixedEffects(auxFixedEffects);
    }
    Superclass::m_StatisticalModel->ComputeCompleteLogLikelihoodGradient(
        Superclass::m_DataSet, auxPopRER, auxIndRER, popGrad, indGrad, gradSqNorms);

    ScalarType deltaF_cur = m_LogLikelihoodTermsHistory[iter - 1].sum() - m_LogLikelihoodTermsHistory[iter].sum();
    ScalarType deltaF_ref = m_LogLikelihoodTermsHistory[iterRef].sum() - m_LogLikelihoodTermsHistory[iter].sum();

    if (fabs(deltaF_cur) < m_AdaptiveTolerance * fabs(deltaF_ref)) {
      std::cout << "Tolerance threshold met. Stopping the optimization process.\n" << std::endl;
      break;
    }

    /* SERIALIZATION */
    if (def::utils::settings.save_state && !(iter % Superclass::m_SaveEveryNIters)) {
      deformation_state << computation_end_state << iter << freezeDirectionCounter << tau
                        << step << gradSqNorms << m_LogLikelihoodTermsHistory
                        << fixedEffects << auxFixedEffects
                        << Superclass::m_PopulationRER << auxPopRER
                        << Superclass::m_IndividualRER << auxIndRER
                        << popGrad << indGrad;
      deformation_state.save_and_reset(def::utils::settings.output_state_filename);
    }
  } // Main loop.

  Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);

  std::cout << "Write output files ...";
  Superclass::m_StatisticalModel->Write(Superclass::m_DataSet,
                                        Superclass::m_PopulationRER, Superclass::m_IndividualRER);
  std::cout << " done." << std::endl;

  /* SERIALIZATION */
  computation_end_state = true;
  if (def::utils::settings.save_state) {
    deformation_state << computation_end_state << iter << freezeDirectionCounter << tau
                      << step << gradSqNorms << m_LogLikelihoodTermsHistory
                      << fixedEffects << auxFixedEffects
                      << Superclass::m_PopulationRER << auxPopRER
                      << Superclass::m_IndividualRER << auxIndRER
                      << popGrad << indGrad;
    deformation_state.save_and_reset(def::utils::settings.output_state_filename);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
FastGradientAscent<ScalarType, Dimension>
::GradientAscentStep(const LinearVariableMapType &auxFixedEffects,
                     const LinearVariableMapType &auxPopRER,
                     const LinearVariablesMapType &auxIndRER,
                     const LinearVariableMapType &popGrad,
                     const LinearVariablesMapType &indGrad,
                     const VectorType &step,
                     LinearVariableMapType &newFixedEffects,
                     LinearVariableMapType &newPopRER,
                     LinearVariablesMapType &newIndRER) {
  newFixedEffects = auxFixedEffects;
  newPopRER = auxPopRER;
  newIndRER = auxIndRER;

  unsigned int k = 0;

  for (auto it = popGrad.begin(); it != popGrad.end(); ++it, ++k) {
    if (newFixedEffects.count(it->first)) {
      newFixedEffects[it->first] = auxFixedEffects.at(it->first) + it->second * step(k);
    } else {
      newPopRER[it->first] = auxPopRER.at(it->first) + it->second * step(k);
    }
  }

  for (auto it = indGrad.begin(); it != indGrad.end(); ++it, ++k) {
    newIndRER[it->first] = auxIndRER.at(it->first) + it->second * step(k);
  }
}

template
class FastGradientAscent<double, 2>;
template
class FastGradientAscent<double, 3>;

