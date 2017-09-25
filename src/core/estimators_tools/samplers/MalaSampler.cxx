/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "MalaSampler.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
MalaSampler<ScalarType, Dimension>
::MalaSampler() : m_Threshold(1000) {
  this->SetMalaType();
}

template<class ScalarType, unsigned int Dimension>
MalaSampler<ScalarType, Dimension>
::~MalaSampler() {}

template<class ScalarType, unsigned int Dimension>
MalaSampler<ScalarType, Dimension>
::MalaSampler(const MalaSampler &other) {}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
MalaSampler<ScalarType, Dimension>
::Sample(LinearVariableMapType &popRER,
         LinearVariablesMapType &indRER,
         VectorType &acceptanceRates) {
  /// Computes the current model and prior terms.
  ProbabilityDistributionMapType popRE, indRE; // "RE" Random Effects.
  Superclass::m_StatisticalModel->GetPopulationRandomEffects(popRE);
  Superclass::m_StatisticalModel->GetIndividualRandomEffects(indRE);

  ScalarType currModelTerm = Superclass::m_StatisticalModel
      ->ComputeModelLogLikelihood(Superclass::m_DataSet, popRER, indRER);

  ScalarType currPriorTerm = 0.0;
  for (auto it = popRER.begin(); it != popRER.end(); ++it) {
    std::shared_ptr<ProbabilityDistributionType> priorRED = popRE[it->first];
    currPriorTerm += priorRED->ComputeLogLikelihood(it->second);
  }
  for (auto it = indRER.begin(); it != indRER.end(); ++it) {
    std::shared_ptr<ProbabilityDistributionType> priorRED = indRE[it->first];
    for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i)
      currPriorTerm += priorRED->ComputeLogLikelihood(it->second[i]);
  }

  /// Draws the candidate and computes the associated model and prior terms.
  VectorType candRER = m_CurrentProposalRED.Sample();

  LinearVariableMapType candPopRER;
  LinearVariablesMapType candIndRER;

  ScalarType candPriorTerm = 0.0;
  unsigned int RECount = 0, scalarCount = 0;
  for (auto it = popRER.begin(); it != popRER.end(); ++it, ++RECount) {
    const unsigned int n = it->second.n_elem();
    candPopRER[it->first] = unvectorize(candRER.subvec(scalarCount, scalarCount + n - 1), m_SizeParametersRER[RECount]);
    scalarCount += n;

    std::shared_ptr<ProbabilityDistributionType> priorRED = popRE[it->first];
    candPriorTerm += priorRED->ComputeLogLikelihood(candPopRER[it->first]);
  }
  for (auto it = indRER.begin(); it != indRER.end(); ++it, ++RECount) {
    LinearVariablesType aux(Superclass::m_NumberOfSubjects);
    const unsigned int n = it->second[0].n_elem();
    std::shared_ptr<ProbabilityDistributionType> priorRED = indRE[it->first];
    for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i) {
      aux[i] = unvectorize(candRER.subvec(scalarCount, scalarCount + n - 1), m_SizeParametersRER[RECount]);
      scalarCount += n;

      candPriorTerm += indRE[it->first]->ComputeLogLikelihood(aux[i]);
    }
    candIndRER[it->first] = aux;
  }

  ScalarType candModelTerm = Superclass::m_StatisticalModel
      ->ComputeModelLogLikelihood(Superclass::m_DataSet, candPopRER, candIndRER);

  /// Computes the candidate proposal distribution and the transition terms.
  VectorType candidateTotalRER;
  NormalDistributionType candidateProposalRED;
  ComputeMalaProposalDistribution(candPopRER, candIndRER, candidateTotalRER, candidateProposalRED);

  ScalarType directTransitionTerm = m_CurrentProposalRED.ComputeLogLikelihood(candidateTotalRER);
  ScalarType inverseTransitionTerm = candidateProposalRED.ComputeLogLikelihood(m_CurrentTotalRER);

  /// Accept-reject.
  acceptanceRates.set_size(popRER.size() + indRER.size());
  UniformDistributionType u;

  ScalarType tau = candPriorTerm + candModelTerm + inverseTransitionTerm
      - currPriorTerm - currModelTerm - directTransitionTerm;

  /// For bug-tracking.
  std::cout << "tau = " << tau
            << "\t[MT = " << candModelTerm - currModelTerm
            << "\tPT = " << candPriorTerm - currPriorTerm
            << "\tTT = " << inverseTransitionTerm - directTransitionTerm
            << "]" << std::endl;

  if (std::log(u.Sample()(0)) > tau) // Reject.
  {
    acceptanceRates.fill(0.0);
  } else                               // Accept.
  {
    popRER = candPopRER;
    indRER = candIndRER;
    m_CurrentTotalRER = candidateTotalRER;
    m_CurrentProposalRED = candidateProposalRED;
    acceptanceRates.fill(100.0);
  }
}

template<class ScalarType, unsigned int Dimension>
void
MalaSampler<ScalarType, Dimension>
::Initialize(LinearVariableMapType const &popRER,
             LinearVariablesMapType const &indRER) {
  /// Initializes the needed size attributes of the random effects realizations.
  m_SizeParametersRER.resize(popRER.size() + indRER.size());
  m_TotalSizeRER = 0;
  unsigned int countRE = 0;
  for (auto it = popRER.begin(); it != popRER.end(); ++it, ++countRE) {
    VectorType aux = vectorize(it->second, m_SizeParametersRER[countRE]);
    m_TotalSizeRER += it->second.n_elem();
  }
  for (auto it = indRER.begin(); it != indRER.end(); ++it) {
    VectorType aux = vectorize(it->second[0], m_SizeParametersRER[countRE]);
    m_TotalSizeRER += it->second[0].n_elem() * Superclass::m_NumberOfSubjects;
  }

  /// Initializes the memory of the current proposal distribution.
  ComputeMalaProposalDistribution(popRER, indRER, m_CurrentTotalRER, m_CurrentProposalRED);
}

template<class ScalarType, unsigned int Dimension>
void
MalaSampler<ScalarType, Dimension>
::AdaptProposalDistributions(const VectorType &detectedAcceptanceRates,
                             const unsigned int &iterationNumber,
                             const bool verbose) {
  const ScalarType goal = Superclass::m_AcceptanceRatesTarget;
  const ScalarType ar = detectedAcceptanceRates(0);

  if (verbose) {
    std::cout << "Scale parameters re-evaluated from ";
    for (auto it = m_Scales.begin(); it != m_Scales.end(); ++it) {
      std::cout << "\t" << it->second << " [" << it->first << "]";
    }
    std::cout << " ..." << std::endl;
  }

  if (ar > goal) {
    for (auto it = m_Scales.begin(); it != m_Scales.end(); ++it) {
      it->second *= 1 + (ar - goal) / ((100 - goal) * sqrt(iterationNumber));
    }
  } else {
    for (auto it = m_Scales.begin(); it != m_Scales.end(); ++it) {
      it->second *= 1 - (goal - ar) / (goal * sqrt(iterationNumber));
    }
  }

  if (verbose) {
    std::cout << "... to ";
    for (auto it = m_Scales.begin(); it != m_Scales.end(); ++it) {
      std::cout << "\t" << it->second << " [" << it->first << "]";
    }
    std::cout << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s)
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
MalaSampler<ScalarType, Dimension>
::ComputeMalaProposalDistribution(LinearVariableMapType const &popRER,
                                  LinearVariablesMapType const &indRER,
                                  VectorType &totRER,
                                  NormalDistributionType &proposalRED) const {
  /// Computes the log-likelihood gradient wrt the random effects.
  LinearVariableMapType popREGrad;
  LinearVariablesMapType indREGrad;
  Superclass::m_StatisticalModel
      ->ComputeLogLikelihoodGradientWrtRandomEffects(Superclass::m_DataSet, popRER, indRER,
                                                     popREGrad, indREGrad);

  /// Vectorizes the random effects realization and the computed gradient.
  totRER.set_size(m_TotalSizeRER);
  totRER.fill(0.0);
  VectorType proposalMean(m_TotalSizeRER, 0.0);
  MatrixType proposalCovarianceSqrt(m_TotalSizeRER, m_TotalSizeRER, 0.0);
  unsigned int RECount = 0, scalarCount = 0;
  for (auto it = popRER.begin(); it != popRER.end(); ++it, ++RECount) {
    totRER.update(vectorize(it->second), scalarCount);

    VectorType drift_block = vectorize(popREGrad[it->first]);
    for (unsigned int k = 0; k < drift_block.size(); ++k)
      drift_block(k) *= m_Threshold / std::max(m_Threshold, std::fabs(drift_block(k)));

    const unsigned int n = it->second.n_elem();
    proposalMean.update(m_Scales.at(it->first) * drift_block, scalarCount);
    proposalCovarianceSqrt.update(diagonal_matrix<ScalarType>(n, sqrt(2 * m_Scales.at(it->first))),
                                  scalarCount, scalarCount);
    scalarCount += n;
  }
  for (auto it = indRER.begin(); it != indRER.end(); ++it, ++RECount) {
    const unsigned int n = it->second[0].n_elem();
    const ScalarType sigma = sqrt(2 * m_Scales.at(it->first));
    for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i) {
      totRER.update(vectorize(it->second[i]), scalarCount);

      VectorType drift_block = vectorize(indREGrad[it->first][i]);
      for (unsigned int k = 0; k < drift_block.size(); ++k)
        drift_block(k) *= m_Threshold / std::max(m_Threshold, std::fabs(drift_block(k)));

      proposalMean.update(m_Scales.at(it->first) * drift_block, scalarCount);
      proposalCovarianceSqrt.update(diagonal_matrix<ScalarType>(n, sigma), scalarCount, scalarCount);
      scalarCount += n;
    }
  }

  /// Sets up the corresponding proposal distribution.
  proposalMean += totRER;

  proposalRED.SetMean(proposalMean);
  proposalRED.SetCovarianceSqrt(proposalCovarianceSqrt);
}

template
class MalaSampler<double, 2>;
template
class MalaSampler<double, 3>;

//#endif /* _MalaSampler_txx */