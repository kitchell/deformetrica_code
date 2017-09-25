/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "GradientAscent.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
GradientAscent<ScalarType, Dimension>
::GradientAscent() : Superclass(), m_InitialStepSize(0.001), m_MultivariateLineSearchFlag(true) {
  this->SetGradientAscentType();
}

template<class ScalarType, unsigned int Dimension>
GradientAscent<ScalarType, Dimension>
::~GradientAscent() {}

template<class ScalarType, unsigned int Dimension>
GradientAscent<ScalarType, Dimension>
::GradientAscent(const GradientAscent &other) {
  m_LogLikelihoodTermsHistory = other.m_LogLikelihoodTermsHistory;
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
GradientAscent<ScalarType, Dimension>
::Update() {

  /// Declare variables.
  LinearVariableMapType fixedEffects, newFixedEffects, newPopRER, popGrad;
  LinearVariablesMapType newIndRER, indGrad;

  unsigned int iter(1), iterRef(0), k(0);
  unsigned long nbParams;
  ScalarType lsqRef;
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

    if (computation_end_state) {
      std::cout << std::endl << "WARNING: Computation was already completed using the current serialize file"
                << std::endl;
    }

    deformation_state >> iter >> step >> m_LogLikelihoodTermsHistory
                      >> fixedEffects >> Superclass::m_PopulationRER >> Superclass::m_IndividualRER
                      >> popGrad >> indGrad;

    Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);
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

    m_LogLikelihoodTermsHistory.resize(Superclass::m_MaxIterations + 1);
    Superclass::m_StatisticalModel->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
        this->m_DataSet, this->m_PopulationRER, this->m_IndividualRER, m_LogLikelihoodTermsHistory[0]);
    lsqRef = m_LogLikelihoodTermsHistory[0].sum();

    Print();

    Superclass::m_StatisticalModel->ComputeCompleteLogLikelihoodGradient(
        Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER, popGrad, indGrad);

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

      GradientAscentStep(fixedEffects, popGrad, indGrad, step, newFixedEffects, newPopRER, newIndRER);
      Superclass::m_StatisticalModel->SetFixedEffects(newFixedEffects);

      bool oob = Superclass::m_StatisticalModel->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
          Superclass::m_DataSet, newPopRER, newIndRER, newLogLikelihoodTerms);
      ScalarType Q = newLogLikelihoodTerms.sum() - lsqRef;

      if (Q >= 0 && !oob) // Easy case : J < Jcurr.
      {
        foundMin = true;
        break;
      } else if (nbParams > 1 && m_MultivariateLineSearchFlag) {
        // All but one step are decreased. All configurations are tried.
        step *= m_AdaptiveShrink;
        std::vector<LinearVariableMapType> newFixedEffects_prop(nbParams);
        std::vector<LinearVariableMapType> newPopRER_prop(nbParams);
        std::vector<LinearVariablesMapType> newIndRER_prop(nbParams);
        std::vector<VectorType> newLogLikelihoodTerms_prop(nbParams);
        std::vector<bool> oob_prop(nbParams);
        VectorType Q(nbParams);

        for (unsigned int k = 0; k < nbParams; ++k) {
          VectorType localStep = step;
          localStep(k) /= m_AdaptiveShrink;

          GradientAscentStep(fixedEffects, popGrad, indGrad, localStep,
                             newFixedEffects_prop[k], newPopRER_prop[k], newIndRER_prop[k]);
          Superclass::m_StatisticalModel->SetFixedEffects(newFixedEffects_prop[k]);

          oob_prop[k] = Superclass::m_StatisticalModel->UpdateFixedEffectsAndComputeCompleteLogLikelihood(
              Superclass::m_DataSet, newPopRER_prop[k], newIndRER_prop[k], newLogLikelihoodTerms_prop[k]);
          Q[k] = newLogLikelihoodTerms_prop[k].sum() - lsqRef;
        }

        unsigned int index = Q.index_max();

        if (Q[index] >= 0 && !oob_prop[index]) {
          newLogLikelihoodTerms = newLogLikelihoodTerms_prop[index];
          newFixedEffects.clear();
          newFixedEffects = newFixedEffects_prop[index];
          newPopRER = newPopRER_prop[index];
          newIndRER = newIndRER_prop[index];
          step(index) /= m_AdaptiveShrink;
          foundMin = true;
          break;
        }
      } else { step *= m_AdaptiveShrink; }
    } // Line search loop.

    if (!foundMin) // Line search loop terminated without finding larger log-likelihood.
    {
      Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);
      std::cout << "Number of line search loops exceeded." << std::endl;
      break;
    }

    m_LogLikelihoodTermsHistory[iter] = newLogLikelihoodTerms;
    fixedEffects = newFixedEffects;
    Superclass::m_PopulationRER = newPopRER;
    Superclass::m_IndividualRER = newIndRER;
    step *= m_AdaptiveExpand;
    lsqRef = newLogLikelihoodTerms.sum();

    /// Displays information about the current state of the algorithm.
    if (!(iter % Superclass::m_PrintEveryNIters)) { Print(); }
    if (!(iter % Superclass::m_SaveEveryNIters)) {
      Superclass::m_StatisticalModel->Write(Superclass::m_DataSet,
                                            Superclass::m_PopulationRER, Superclass::m_IndividualRER);
    }

    ScalarType deltaF_cur = m_LogLikelihoodTermsHistory[iter - 1].sum() - m_LogLikelihoodTermsHistory[iter].sum();
    ScalarType deltaF_ref = m_LogLikelihoodTermsHistory[iterRef].sum() - m_LogLikelihoodTermsHistory[iter].sum();

    if (fabs(deltaF_cur) < m_AdaptiveTolerance * fabs(deltaF_ref)) {
      std::cout << "Tolerance threshold met. Stopping the optimization process.\n" << std::endl;
      break;
    }

    Superclass::m_StatisticalModel->SetFixedEffects(fixedEffects);
    Superclass::m_StatisticalModel->ComputeCompleteLogLikelihoodGradient(
        Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER, popGrad, indGrad);

    /* SERIALIZATION */
    if (def::utils::settings.save_state && !(iter % Superclass::m_SaveEveryNIters)) {
      deformation_state << computation_end_state << iter << step << m_LogLikelihoodTermsHistory
                        << fixedEffects << Superclass::m_PopulationRER << Superclass::m_IndividualRER
                        << popGrad << indGrad;
      deformation_state.save_and_reset(def::utils::settings.output_state_filename);
    }
  }

  std::cout << "Write output files ...";
  Superclass::m_StatisticalModel->Write(
      Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER);
  std::cout << " done." << std::endl;

  /* SERIALIZATION */
  computation_end_state = true;
  if (def::utils::settings.save_state) {
    deformation_state << computation_end_state << iter << step << m_LogLikelihoodTermsHistory
                      << fixedEffects << Superclass::m_PopulationRER << Superclass::m_IndividualRER
                      << popGrad << indGrad;
    deformation_state.save_and_reset(def::utils::settings.output_state_filename);
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
GradientAscent<ScalarType, Dimension>
::GradientAscentStep(const LinearVariableMapType &fixedEffects,
                     const LinearVariableMapType &popGrad,
                     const LinearVariablesMapType &indGrad,
                     const VectorType &step,
                     LinearVariableMapType &newFixedEffects,
                     LinearVariableMapType &newPopRER,
                     LinearVariablesMapType &newIndRER) {
  newFixedEffects.clear();
  newFixedEffects = fixedEffects;
  newPopRER = Superclass::m_PopulationRER;
  newIndRER = Superclass::m_IndividualRER;

  unsigned int k = 0;

  for (auto it = popGrad.begin(); it != popGrad.end(); ++it, ++k) {
    if (newFixedEffects.count(it->first)) {
      newFixedEffects[it->first] = fixedEffects.at(it->first) + it->second * step(k);
    } else {
      newPopRER[it->first] = Superclass::m_PopulationRER[it->first] + it->second * step(k);
    }
  }

  for (auto it = indGrad.begin(); it != indGrad.end(); ++it, ++k) {
    newIndRER[it->first] = Superclass::m_IndividualRER[it->first] + it->second * step(k);
  }
}

template
class GradientAscent<double, 2>;
template
class GradientAscent<double, 3>;

