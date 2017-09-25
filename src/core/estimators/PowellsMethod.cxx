/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "PowellsMethod.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
PowellsMethod<ScalarType, Dimension>
::PowellsMethod() : Superclass(), m_MaxLineSearchIterations(10) {
  this->SetPowellsMethodType();
}

template<class ScalarType, unsigned int Dimension>
PowellsMethod<ScalarType, Dimension>
::~PowellsMethod() {}

template<class ScalarType, unsigned int Dimension>
PowellsMethod<ScalarType, Dimension>
::PowellsMethod(const PowellsMethod &other) {
  m_MaxLineSearchIterations = other.m_MaxLineSearchIterations;
  m_InitialUncertainty = other.m_InitialUncertainty;
  m_Structure = other.m_Structure;
  m_Keys = other.m_Keys;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
PowellsMethod<ScalarType, Dimension>
::Update() {

  /// Initialization.
  LinearVariableMapType fixedEffects;
  Superclass::m_StatisticalModel->GetFixedEffects(fixedEffects);

  VectorType initialPosition = fixedEffects.vectorize(m_Structure, m_Keys);
  VectorType currentPosition = initialPosition;
  const unsigned int problemDimension = currentPosition.size();

  VectorType functionalValues(problemDimension + 1, 0.0);
  VectorType increments(problemDimension, 0.0);
  functionalValues(0) = Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(
      Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER);

  VectorType initialUncertainty = LinearizeInitialUncertainty(fixedEffects);
  std::vector<VectorType> currentDirections(problemDimension);
  for (unsigned int k = 0; k < problemDimension; ++k) {
    currentDirections[k].set_size(problemDimension);
    currentDirections[k].fill(0.0);
    currentDirections[k](k) = initialUncertainty(k);
  }

  /// Main loop.
  for (unsigned int iter = 1; iter < Superclass::m_MaxIterations + 1; ++iter) {
    /// Memory of the initial position.
    initialPosition = currentPosition;

    /// Iteratively move along the current set of directions.
    for (unsigned int k = 0; k < problemDimension; ++k) {
      functionalValues(k + 1) = LineSearchStep(currentDirections[k], functionalValues(k), currentPosition);
      increments(k) = functionalValues(k + 1) - functionalValues(k);
    }

    /// Direction of maximum improvement.
    const unsigned int bestDirectionIndex = increments.index_max();
    const ScalarType maxIncrement = increments(bestDirectionIndex);

    /// Functional value in the extended direction.
    const VectorType extendedDirection = 2 * currentPosition - initialPosition;
    const ScalarType extendedDirectionFunctionalValue = ComputeFunctional(extendedDirection);

    /// Check if a last line minimization in the main direction is interesting before looping.
    const ScalarType f0 = functionalValues(0);
    const ScalarType fn = functionalValues(problemDimension);
    const ScalarType fe = extendedDirectionFunctionalValue;
    if (fe < f0 || 2 * (f0 - 2 * fn + fe) * pow(f0 - fn - maxIncrement, 2) < maxIncrement * pow(f0 - fe, 2)) {
      Superclass::m_StatisticalModel->SetFixedEffects(unvectorize(currentPosition, m_Structure, m_Keys));
      functionalValues(0) = fn;

    } else {
      currentDirections[bestDirectionIndex] = currentPosition - initialPosition;
      functionalValues(0) = LineSearchStep(currentDirections[bestDirectionIndex], fn, currentPosition);
      Superclass::m_StatisticalModel->SetFixedEffects(unvectorize(currentPosition, m_Structure, m_Keys));
    }
  }

  /// Write the outputs.
  std::cout << "Write output files ...";
  Superclass::m_StatisticalModel->Write(
      Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER);
  std::cout << " done." << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
ScalarType
PowellsMethod<ScalarType, Dimension>
::ComputeFunctional(const VectorType &position) {
  Superclass::m_StatisticalModel->SetFixedEffects(unvectorize(position, m_Structure, m_Keys));

  return Superclass::m_StatisticalModel->ComputeCompleteLogLikelihood(
      Superclass::m_DataSet, Superclass::m_PopulationRER, Superclass::m_IndividualRER);
}

template<class ScalarType, unsigned int Dimension>
ScalarType
PowellsMethod<ScalarType, Dimension>
::LineSearchStep(const VectorType &direction,
                 const ScalarType functional,
                 VectorType &position) {

  const ScalarType searchRatio = 0.3819660115;
  const ScalarType targetAccuracy = 1e-3;

  /// Initialization.
  ScalarType leftFrontier = - 1.0;
  ScalarType rightFrontier = 1.0;
  ScalarType segmentLength = 2.0;

  ScalarType leftFrontierValue = ComputeFunctional(position + leftFrontier * direction);
  ScalarType rightFrontierValue = ComputeFunctional(position + rightFrontier * direction);

  /// First iteration, taking advantage of the known functional value at the current position.
  ScalarType leftInterior = 0.0;
  ScalarType rightInterior = 0.5;

  ScalarType leftInteriorValue = functional;
  ScalarType rightInteriorValue = ComputeFunctional(position + rightInterior * direction);

  if (leftInteriorValue < rightInteriorValue) {
    rightFrontier = rightInterior;
    segmentLength = rightFrontier - leftFrontier;

    rightFrontierValue = rightInteriorValue;

    leftInterior = leftFrontier + searchRatio * segmentLength;
    rightInterior = rightFrontier - searchRatio * segmentLength;

    ScalarType leftInteriorValue = ComputeFunctional(position + leftInterior * direction);
    ScalarType rightInteriorValue = ComputeFunctional(position + rightInterior * direction);

  } else {
    leftFrontier = leftInterior;
    segmentLength = rightFrontier - leftFrontier;

    leftFrontierValue = leftInteriorValue;

    leftInterior = leftFrontier + searchRatio * segmentLength;
    rightInterior = rightFrontier - searchRatio * segmentLength;

    ScalarType leftInteriorValue = ComputeFunctional(position + leftInterior * direction);
    ScalarType rightInteriorValue = ComputeFunctional(position + rightInterior * direction);
  }

  /// Main loop.
  for (unsigned int iter = 0; iter < m_MaxLineSearchIterations; ++iter) {
    if (leftInteriorValue < rightInteriorValue) {
      rightFrontier = rightInterior;
      segmentLength = rightFrontier - leftFrontier;
      if (segmentLength < targetAccuracy) { break; }

      rightFrontierValue = rightInteriorValue;

      leftInterior = leftFrontier + searchRatio * segmentLength;
      rightInterior = leftInterior;

      ScalarType leftInteriorValue = ComputeFunctional(position + leftInterior * direction);
      ScalarType rightInteriorValue = leftInteriorValue;

    } else {
      leftFrontier = leftInterior;
      segmentLength = rightFrontier - leftFrontier;
      if (segmentLength < targetAccuracy) { break; }

      leftFrontierValue = leftInteriorValue;

      leftInterior = rightInterior;
      rightInterior = rightFrontier - searchRatio * segmentLength;

      ScalarType leftInteriorValue = ComputeFunctional(position + leftInterior * direction);
      ScalarType rightInteriorValue = ComputeFunctional(position + rightInterior * direction);
    }
  }

  /// Return result.
  position += 0.5 * (leftFrontier + rightFrontier) * direction;
  return ComputeFunctional(position);
}

template<class ScalarType, unsigned int Dimension>
typename def::algebra::VectorType
PowellsMethod<ScalarType, Dimension>
::LinearizeInitialUncertainty(const LinearVariableMapType &modelFixedEffects) const {

  VectorType out;
  for (auto it = modelFixedEffects.begin(); it != modelFixedEffects.end(); ++it) {
    out.push_back(recast<ScalarType>(m_InitialUncertainty.at(it->first)) * VectorType(it->second.n_elem(), 1.0));
  }

  return out;
}

template
class PowellsMethod<double, 2>;
template
class PowellsMethod<double, 3>;

