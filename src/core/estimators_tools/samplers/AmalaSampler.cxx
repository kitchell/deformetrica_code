/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/



#include "AmalaSampler.h"

/// For bug-tracking.
#include "MatrixDLM.h"

/// For bug-tracking.
template<class ScalarType, unsigned int Dimension>
unsigned int AmalaSampler<ScalarType, Dimension>::m_Count = 0;


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
AmalaSampler<ScalarType, Dimension>
::AmalaSampler() : m_Threshold(1000) {
  this->SetAmalaType();
}

template<class ScalarType, unsigned int Dimension>
AmalaSampler<ScalarType, Dimension>
::~AmalaSampler() {}

template<class ScalarType, unsigned int Dimension>
AmalaSampler<ScalarType, Dimension>
::AmalaSampler(const AmalaSampler &other) {}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
AmalaSampler<ScalarType, Dimension>
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

//    /// For bug-tracking.
//    std::cout << "Xk_250 = " << m_CurrentTotalRER(250) << " ; Xc_250 = " << candRER(250) << std::endl;

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
  ComputeAmalaProposalDistribution(candPopRER, candIndRER, candidateTotalRER, candidateProposalRED);

  ScalarType directTransitionTerm = m_CurrentProposalRED.ComputeLogLikelihood(candidateTotalRER);
  ScalarType inverseTransitionTerm = candidateProposalRED.ComputeLogLikelihood(m_CurrentTotalRER);

//    /// For bug-tracking.
//    ScalarType g = m_CurrentProposalRED.ComputeLogLikelihood(m_CurrentTotalRER);
//    ScalarType h = candidateProposalRED.ComputeLogLikelihood(candidateTotalRER);
//
//    VectorType mCurr = m_CurrentProposalRED.GetMean();
//    VectorType mCand = candidateProposalRED.GetMean();
//    MatrixType cCurr = m_CurrentProposalRED.GetCovariance();
//    MatrixType ciCurr = m_CurrentProposalRED.GetCovarianceInverse();
//    MatrixType cCand = candidateProposalRED.GetCovariance();
//    MatrixType ciCand = candidateProposalRED.GetCovarianceInverse();
//
//    VectorType a  = mCand - mCurr;
//    VectorType aa = mCurr - mCand;
//    VectorType b  = ciCurr * a;
//    VectorType bb = ciCand * aa;
//    ScalarType c  = dot_product(a, b);
//    ScalarType cc = dot_product(aa, bb);
//    ScalarType d  = log_det(cCurr);
//    ScalarType dd = log_det(cCand);
//
//    ScalarType e = det(cCurr);
//    ScalarType f = det(ciCurr);

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
//        /// For bug-tracking.
//        std::ostringstream oss5;
//        oss5 << "0_Momenta1_before.txt" << std::ends;
//        writeMatrixDLM<ScalarType>(oss5.str().c_str(), recast<MatrixType>(indRER["Momenta"][0]));

    popRER = candPopRER;
    indRER = candIndRER;
    m_CurrentTotalRER = candidateTotalRER;
    m_CurrentProposalRED = candidateProposalRED;
    acceptanceRates.fill(100.0);

//        /// For bug-tracking.
//        std::ostringstream oss6;
//        oss6 << "0_Momenta1_after.txt" << std::ends;
//        writeMatrixDLM<ScalarType>(oss6.str().c_str(), recast<MatrixType>(indRER["Momenta"][0]));

//        /// For bug-tracking.
//        std::cout << "Variation of complete log-likelihood = "
//                  << candPriorTerm + candModelTerm - currPriorTerm - currModelTerm << std::endl;
  }
}

template<class ScalarType, unsigned int Dimension>
void
AmalaSampler<ScalarType, Dimension>
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
  ComputeAmalaProposalDistribution(popRER, indRER, m_CurrentTotalRER, m_CurrentProposalRED);
}

template<class ScalarType, unsigned int Dimension>
void
AmalaSampler<ScalarType, Dimension>
::AdaptProposalDistributions(const VectorType &detectedAcceptanceRates,
                             const unsigned int &iterationNumber,
                             const bool verbose) {
//    const ScalarType goal = Superclass::m_AcceptanceRatesTarget;
//    const ScalarType ar = detectedAcceptanceRates(0);
//
//    std::cout << "Scale parameter re-evaluated from " << m_Scale;
//
//    if (ar > goal)
//        m_Scale *= 1 + (ar - goal) / ((100 - goal) * sqrt(iterationNumber));
//    else
//        m_Scale *= 1 - (goal - ar) / (goal * sqrt(iterationNumber));
//
//    std::cout << " to " << m_Scale << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s)
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
AmalaSampler<ScalarType, Dimension>
::ComputeAmalaProposalDistribution(LinearVariableMapType const &popRER,
                                   LinearVariablesMapType const &indRER,
                                   VectorType &totRER,
                                   NormalDistributionType &proposalRED) const {
  /// Computes the log-likelihood gradient wrt the random effects.
  LinearVariableMapType popREGrad;
  LinearVariablesMapType indREGrad;
  Superclass::m_StatisticalModel
      ->ComputeCompleteLogLikelihoodGradient(Superclass::m_DataSet, popRER, indRER,
                                             popREGrad, indREGrad);

  /// Vectorizes the random effects realization and the computed gradient.
  totRER.set_size(m_TotalSizeRER);
  totRER.fill(0.0);
  VectorType proposalMean(m_TotalSizeRER, 0.0);
  MatrixType proposalCovariance(m_TotalSizeRER, m_TotalSizeRER, 0.0);
  unsigned int RECount = 0, scalarCount = 0;
  for (auto it = popRER.begin(); it != popRER.end(); ++it, ++RECount) {
    totRER.update(vectorize(it->second), scalarCount);

    MatrixType drift_block = vectorize(popREGrad[it->first]);
    for (unsigned int k = 0; k < drift_block.size(); ++k)
      drift_block(k, 0) *= m_Threshold / std::max(m_Threshold, std::fabs(drift_block(k, 0)));

    const unsigned int n = it->second.n_elem();
    proposalMean.update(m_MeanScales.at(it->first) * drift_block.get_column(0), scalarCount);
    proposalCovariance.update(m_CovarianceScales.at(it->first)
                                  * (diagonal_matrix<ScalarType>(n, m_Regularizations.at(it->first))
                                      + drift_block * drift_block.transpose()), scalarCount, scalarCount);
    scalarCount += n;
  }
  for (auto it = indRER.begin(); it != indRER.end(); ++it, ++RECount) {
    const unsigned int n = it->second[0].n_elem();
    for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i) {
      totRER.update(vectorize(it->second[i]), scalarCount);

      MatrixType drift_block = vectorize(indREGrad[it->first][i]);
      for (unsigned int k = 0; k < drift_block.size(); ++k)
        drift_block(k, 0) *= m_Threshold / std::max(m_Threshold, std::fabs(drift_block(k, 0)));

      proposalMean.update(m_MeanScales.at(it->first) * drift_block.get_column(0), scalarCount);
      proposalCovariance.update(m_CovarianceScales.at(it->first)
                                    * (diagonal_matrix<ScalarType>(n, m_Regularizations.at(it->first))
                                        + drift_block * drift_block.transpose()), scalarCount, scalarCount);
      scalarCount += n;
    }
  }

  /// Sets up the corresponding proposal distribution.
  proposalMean += totRER;

//    /// For bug-tracking (full correlation).
//    MatrixType drift = (proposalMean - totRER) / m_Scale;
//    proposalCovariance = drift * drift.transpose();

//    /// For bug-tracking (only diagonal elements).
//    VectorType drift = (proposalMean - totRER) / m_Scale;
//    proposalCovariance.fill(0.0);
//    ++m_Count;
//    ScalarType lambda = std::max(0.5, 1 - 0.02 * m_Count);
//    std::cout << "lambda = " << lambda << std::endl;
//    for (unsigned int k = 0 ; k < m_TotalSizeRER ; ++k)
//        proposalCovariance(k, k) = sqrt(m_Scale * ((1 - lambda) * drift(k) * drift(k)
//                                                   + std::max(m_Regularization,  lambda)));
//    proposalRED.SetCovarianceSqrt(proposalCovariance);

//    proposalCovariance += diagonal_matrix(scalarCount, m_Regularization);
//    proposalCovariance *= m_Scale;

  proposalRED.SetMean(proposalMean);
  proposalRED.SetCovariance(proposalCovariance);
}

template
class AmalaSampler<double, 2>;
template
class AmalaSampler<double, 3>;
