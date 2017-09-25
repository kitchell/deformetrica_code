/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "SrwMhwgSampler.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
SrwMhwgSampler<ScalarType, Dimension>
::SrwMhwgSampler() {
  this->SetSrwMhwgType();
}

template<class ScalarType, unsigned int Dimension>
SrwMhwgSampler<ScalarType, Dimension>
::~SrwMhwgSampler() {}

template<class ScalarType, unsigned int Dimension>
SrwMhwgSampler<ScalarType, Dimension>
::SrwMhwgSampler(const SrwMhwgSampler &other) {}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
SrwMhwgSampler<ScalarType, Dimension>
::Sample(LinearVariableMapType &popRER,
         LinearVariablesMapType &indRER,
         VectorType &acceptanceRates,
         const ScalarType &temperature) {
//  /// For bug-tracking.
//  ScalarType cll_before = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(Superclass::m_DataSet,
//                                                                                       popRER, indRER);
//  ScalarType deltaCLL = 0.0;

  UniformDistributionType u;
  ProbabilityDistributionMapType popRE, indRE; // "RE" Random Effects.
  Superclass::m_StatisticalModel->GetPopulationRandomEffects(popRE);
  Superclass::m_StatisticalModel->GetIndividualRandomEffects(indRE);

  // Initialization of the memory of the current model terms.
  // The contribution of each subject is stored independently.
  std::vector<ScalarType> currModelTerms;
  Superclass::m_StatisticalModel->ComputeModelLogLikelihood(
      Superclass::m_DataSet, popRER, indRER, temperature, "All", currModelTerms);

  /// [std::map loop] Loop on the population variables (typically control points, template).
  unsigned int popRECount = 0;
  for (auto it = popRER.begin(); it != popRER.end(); ++it, ++popRECount) {
//    /// For bug-tracking.
//    std::cout << "\n[" << it->first << "]" << std::endl;

    // RED (Random Effect Distribution) initialization.
    std::shared_ptr<ProbabilityDistributionType> priorRED = popRE[it->first];
    std::shared_ptr<AbstractNormalDistributionType> propoRED = m_PopulationProposalDistribution.at(it->first);

    // RER (Random Effect Realization) initialization.
    std::vector<unsigned int> sizeParams;
    VectorType currRER = vectorize(it->second, sizeParams);

    /// Miscellaneous initializations.
    ScalarType currPriorTerm = priorRED->ComputeLogLikelihood(it->second, temperature);
    std::vector<ScalarType> candModelTerms;
    acceptanceRates(popRECount) = 0.0;

    /// [VectorType loop] Exhaustive loop on every scalar value of the current and candidate RER.
    unsigned int blockCount = 0;
    for (unsigned int l = 0; l < currRER.size(); ++blockCount) {
      /// Draw the candidate (Cand).
      unsigned int blockSize = propoRED->GetMean().size();
      VectorType currRER_block_memory(blockSize, 0.0);
      for (unsigned int m = 0; m < blockSize;) {
        currRER_block_memory(m) = currRER(l + m);

        ++m;
        if (l + m == currRER.size()) {
          currRER_block_memory.set_size(m);
          blockSize = m;
          break;
        }
      }
      propoRED->SetMean(currRER_block_memory);
      VectorType candRER_block = propoRED->Sample();
      candRER_block.set_size(blockSize);

      /// Evaluate the current part (Curr).
      ScalarType currModelTerm = std::accumulate(currModelTerms.begin(), currModelTerms.end(), 0.0);

      /// Evaluate the candidate part (Cand).
      for (unsigned int m = 0; m < blockSize; ++m)
        currRER(l + m) = candRER_block(m); // Accept

//      /// For bug-tracking.
//      ScalarType before = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(
//          Superclass::m_DataSet, popRER, indRER);

      it->second = unvectorize(currRER, sizeParams);
      ScalarType candPriorTerm = priorRED->ComputeLogLikelihood(currRER, temperature);
      ScalarType candModelTerm = Superclass::m_StatisticalModel->ComputeModelLogLikelihood(
          Superclass::m_DataSet, popRER, indRER, temperature, it->first, candModelTerms);

      /// Accept-reject.
      ScalarType tau = candPriorTerm + candModelTerm - currPriorTerm - currModelTerm;
      if (std::log(u.Sample()(0)) > tau)     // Reject.
      {
        for (unsigned int m = 0; m < blockSize; ++m) { currRER(l + m) = currRER_block_memory(m); }
        it->second = unvectorize(currRER, sizeParams);
        Superclass::m_StatisticalModel->RecoverMemorizedState();
      } else {
        /// For bug-tracking.
//        deltaCLL += candPriorTerm + candModelTerm - currPriorTerm - currModelTerm;
//
//        ScalarType after = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(
//            Superclass::m_DataSet, popRER, indRER);
//        std::cout << "We should have delta model term = " << candModelTerm - currModelTerm << std::endl;
//        std::cout << "We should have delta prior term = " << candPriorTerm - currPriorTerm << std::endl;
//        std::cout << "We should have " << candPriorTerm + candModelTerm - currPriorTerm - currModelTerm
//                  << " = " << after - before;
//        if (std::fabs(after - before - (candPriorTerm + candModelTerm - currPriorTerm - currModelTerm)) < 1e-5) {
//          std::cout << " >> OK" << std::endl;
//        } else { std::cout << " >> NOT OK" << std::endl; }
//        assert(std::fabs(after - before - (candPriorTerm + candModelTerm - currPriorTerm - currModelTerm)) < 1e-5);

        currModelTerms = candModelTerms;
        currPriorTerm = candPriorTerm;
        acceptanceRates(popRECount) += 1.0;
      }

      /// Move on to the next block.
      l += blockSize;
    }
    acceptanceRates(popRECount) *= 100.0 / blockCount;
  }

  /// [std::map loop] Loop on the individual variables (typically momentas). Parallelization ?? TODO.
  unsigned int indRECount = 0;
  for (auto it = indRER.begin(); it != indRER.end(); ++it, ++indRECount) {
//    /// For bug-tracking.
//    std::cout << "\n[" << it->first << "]" << std::endl;

    acceptanceRates(popRECount + indRECount) = 0.0;

    // RED (Random Effect Distribution).
    std::shared_ptr<ProbabilityDistributionType> priorRED = indRE[it->first];
    std::shared_ptr<AbstractNormalDistributionType> propoRED = m_IndividualProposalDistribution.at(it->first);

    /// [std::vector loop] Loop on the individuals i.
    unsigned int blockCount = 0;
    for (unsigned int l = 0; l < indRER[it->first][0].n_elem(); ++blockCount) {
      unsigned int blockSize = propoRED->GetMean().size();
      std::vector<ScalarType> currPriorTerms(Superclass::m_NumberOfSubjects);
      std::vector<ScalarType> candPriorTerms(Superclass::m_NumberOfSubjects);
      std::vector<ScalarType> candModelTerms(Superclass::m_NumberOfSubjects);
      std::vector<VectorType> currRERs(Superclass::m_NumberOfSubjects);
      std::vector<VectorType> currRERs_block_memory(Superclass::m_NumberOfSubjects);
      std::vector<VectorType> candRERs_block(Superclass::m_NumberOfSubjects);
      std::vector<unsigned int> sizeParams;

      for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i) {
        currPriorTerms[i] = priorRED->ComputeLogLikelihood(it->second[i], temperature);

        // RER (Random Effect Realization) ; Curr (Current).
        currRERs[i] = vectorize(indRER[it->first][i], sizeParams);

        /// [VectorType loop] Exhaustive loop on every scalar value of the current and candidate RER.
        /// Draw the candidate (Cand).
        currRERs_block_memory[i].set_size(blockSize);
        currRERs_block_memory[i].fill(0.0);
        for (unsigned int m = 0; m < blockSize;) {
          currRERs_block_memory[i](m) = currRERs[i](l + m);
          ++m;
          if (l + m == currRERs[i].size()) {
            currRERs_block_memory[i].set_size(m);
            blockSize = m;
            break;
          }
        }
        propoRED->SetMean(currRERs_block_memory[i]);
        VectorType candRER_i_block = propoRED->Sample();
        candRER_i_block.set_size(blockSize);

        /// Evaluate the candidate part (Cand).
        for (unsigned int m = 0; m < blockSize; ++m)
          currRERs[i](l + m) = candRER_i_block(m); // Accept

//        /// For bug-tracking.
//        ScalarType before = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(
//            Superclass::m_DataSet, popRER, indRER);

        it->second[i] = unvectorize(currRERs[i], sizeParams);
        candPriorTerms[i] = priorRED->ComputeLogLikelihood(currRERs[i], temperature);
      }

      ScalarType candModelTerm = Superclass::m_StatisticalModel->ComputeModelLogLikelihood(
          Superclass::m_DataSet, popRER, indRER, temperature, it->first, candModelTerms);

      for (unsigned long i = 0; i < Superclass::m_NumberOfSubjects; ++i) {
        /// Accept-reject.
        ScalarType tau = candPriorTerms[i] + candModelTerms[i] - currPriorTerms[i] - currModelTerms[i];
        if (std::log(u.Sample()(0)) > tau)       // Reject.
        {
          for (unsigned int m = 0; m < blockSize; ++m) { currRERs[i](l + m) = currRERs_block_memory[i](m); }
          it->second[i] = unvectorize(currRERs[i], sizeParams);
          Superclass::m_StatisticalModel->RecoverMemorizedState();
        } else {
//          /// For bug-tracking.
//          deltaCLL += candPriorTerms[i] + candModelTerms[i] - currPriorTerms[i] - currModelTerms[i];
//
//          ScalarType after = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(
//              Superclass::m_DataSet, popRER, indRER);
//          std::cout << "We should have delta model term = " << candModelTerms[i] - currModelTerms[i] << std::endl;
//          std::cout << "We should have delta prior term = " << candPriorTerms[i] - currPriorTerms[i] << std::endl;
//          std::cout << "We should have " << candPriorTerms[i] + candModelTerms[i] - currPriorTerms[i] - currModelTerms[i]
//                    << " = " << after - before;
//          if (std::fabs(after - before - (candPriorTerms[i] + candModelTerms[i] - currPriorTerms[i] - currModelTerms[i])) < 1e-5) {
//            std::cout << " >> OK \n" << std::endl;
//          } else { std::cout << " >> NOT OK \n" << std::endl; }
//          assert(std::fabs(after - before - (candPriorTerms[i] + candModelTerms[i] - currPriorTerms[i] - currModelTerms[i])) < 1e-5);

          currModelTerms[i] = candModelTerms[i];
          currPriorTerms[i] = candPriorTerms[i];
          acceptanceRates(popRECount + indRECount) += 1.0;
        }
      }
      /// Move on to the next block.
      l += blockSize;
    }
    acceptanceRates(popRECount + indRECount) *= 100.0 / (blockCount * Superclass::m_NumberOfSubjects);
  }

  /// For bug-tracking.
//  std::cout << "Variation of complete log-likelihood = " << deltaCLL << std::endl;
//  ScalarType cll_after = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(Superclass::m_DataSet,
//                                                                                      popRER, indRER);
//  assert(std::fabs(cll_after - cll_before - deltaCLL) < 1e-10);
}

template<class ScalarType, unsigned int Dimension>
void
SrwMhwgSampler<ScalarType, Dimension>
::AdaptProposalDistributions(const VectorType &detectedAcceptanceRates,
                             const unsigned int &iterationNumber,
                             const bool verbose) {
  const ScalarType goal = Superclass::m_AcceptanceRatesTarget;

  if (verbose) { std::cout << ">> Proposal std re-evaluated from :"; }

  unsigned int k = 0;
  for (auto it = m_PopulationProposalDistribution.begin(); it != m_PopulationProposalDistribution.end(); ++it, ++k) {
    const ScalarType ar = detectedAcceptanceRates(k);
    MatrixType covsqrt = it->second->GetCovarianceSqrt();

    if (verbose) {
      std::cout << std::endl;
      std::cout << "\t\t" << covsqrt(0, 0);
    }

    if (ar > goal) { covsqrt *= 1 + (ar - goal) / ((100 - goal) * sqrt(iterationNumber)); }
    else { covsqrt *= 1 - (goal - ar) / (goal * sqrt(iterationNumber)); }

    if (verbose) { std::cout << " to " << covsqrt(0, 0) << " \t[" << it->first << "]"; }

    it->second->SetCovarianceSqrt(covsqrt);
  }
  for (auto it = m_IndividualProposalDistribution.begin(); it != m_IndividualProposalDistribution.end(); ++it, ++k) {
    const ScalarType ar = detectedAcceptanceRates(k);
    MatrixType covsqrt = it->second->GetCovarianceSqrt();

    if (verbose) {
      std::cout << std::endl;
      std::cout << "\t\t" << covsqrt(0, 0);
    }

    if (ar > goal) { covsqrt *= 1 + (ar - goal) / ((100 - goal) * sqrt(iterationNumber)); }
    else { covsqrt *= 1 - (goal - ar) / (goal * sqrt(iterationNumber)); }

    if (verbose) { std::cout << " to " << covsqrt(0, 0) << " \t[" << it->first << "]"; }

    it->second->SetCovarianceSqrt(covsqrt);
  }

  if (verbose) { std::cout << std::endl; }
}

template<class ScalarType, unsigned int Dimension>
void
SrwMhwgSampler<ScalarType, Dimension>
::SaveState(def::utils::DeformationState &state) const {
  LinearVariableMapType populationProposalCovarianceSqrt, individualProposalCovarianceSqrt;

  for (auto it = m_PopulationProposalDistribution.begin(); it != m_PopulationProposalDistribution.end(); ++it) {
    populationProposalCovarianceSqrt[it->first] = it->second->GetCovarianceSqrt();
  }

  for (auto it = m_IndividualProposalDistribution.begin(); it != m_IndividualProposalDistribution.end(); ++it) {
    individualProposalCovarianceSqrt[it->first] = it->second->GetCovarianceSqrt();
  }

  state << populationProposalCovarianceSqrt << individualProposalCovarianceSqrt;
}

template<class ScalarType, unsigned int Dimension>
void
SrwMhwgSampler<ScalarType, Dimension>
::RecoverState(def::utils::DeformationState &state) {
  LinearVariableMapType populationProposalCovarianceSqrt, individualProposalCovarianceSqrt;
  state >> populationProposalCovarianceSqrt >> individualProposalCovarianceSqrt;

  for (auto it = m_PopulationProposalDistribution.begin(); it != m_PopulationProposalDistribution.end(); ++it) {
    it->second->SetCovarianceSqrt(recast<MatrixType>(populationProposalCovarianceSqrt[it->first]));
  }

  for (auto it = m_IndividualProposalDistribution.begin(); it != m_IndividualProposalDistribution.end(); ++it) {
    it->second->SetCovarianceSqrt(recast<MatrixType>(individualProposalCovarianceSqrt[it->first]));
  }
}

template
class SrwMhwgSampler<double, 2>;
template
class SrwMhwgSampler<double, 3>;

