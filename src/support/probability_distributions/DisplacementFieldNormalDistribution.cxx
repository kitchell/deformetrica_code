/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "DisplacementFieldNormalDistribution.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
DisplacementFieldNormalDistribution<ScalarType>
::DisplacementFieldNormalDistribution()
    : m_KernelType(null), m_KernelWidth(0), m_SubsamplingStepSize(0), m_ShapeDimension(0) {}

template<class ScalarType>
DisplacementFieldNormalDistribution<ScalarType>
::~DisplacementFieldNormalDistribution() {}

template<class ScalarType>
DisplacementFieldNormalDistribution<ScalarType>
::DisplacementFieldNormalDistribution(const DisplacementFieldNormalDistribution &other) : Superclass(other) {
  m_MomentaDistribution = other.m_MomentaDistribution;
  m_KernelType = other.m_KernelType;
  m_KernelWidth = other.m_KernelWidth;
  m_SubsamplingStepSize = other.m_SubsamplingStepSize;
  m_ShapeDimension = other.m_ShapeDimension;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType>
VectorType
DisplacementFieldNormalDistribution<ScalarType>
::Sample() const {
  /// Initialization of the random number generator.
  std::random_device rd;
  std::mt19937_64 generator(rd());
  std::normal_distribution<ScalarType> stdNormal(0, 1);

  /// Further initializations.
  const unsigned int numberOfShapePoints = Superclass::m_Mean.size() / m_ShapeDimension;
  const unsigned int numberOfShapeControlPoints = 1 + (numberOfShapePoints - 1) / m_SubsamplingStepSize;
  const MatrixType shapePoints = Superclass::m_Mean.unvectorize(numberOfShapePoints, m_ShapeDimension);

  /// List the shape control points.
  MatrixType shapeControlPoints(numberOfShapeControlPoints, m_ShapeDimension, 0.0);
  for (unsigned int k = 0; k < numberOfShapeControlPoints; ++k) {
    shapeControlPoints.set_row(k, shapePoints.get_row(k * m_SubsamplingStepSize));
  }

  /// Draw the momenta.
  const MatrixType momenta = m_MomentaDistribution.Sample().unvectorize(numberOfShapeControlPoints, m_ShapeDimension);

  /// Convolve, sum, and return.
  MatrixType displacementField;
  if (m_ShapeDimension == 2) {

    typedef KernelFactory<ScalarType, 2> KernelFactoryType;
    typedef typename KernelFactoryType::KernelBaseType KernelType;
    KernelFactoryType *kernelFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernel = kernelFactory->CreateKernelObject(m_KernelType);
    kernel->SetKernelWidth(m_KernelWidth);
    kernel->SetSources(shapeControlPoints);
    kernel->SetWeights(momenta);
    displacementField = kernel->Convolve(shapePoints);

  } else if (m_ShapeDimension == 3) {

    typedef KernelFactory<ScalarType, 3> KernelFactoryType;
    typedef typename KernelFactoryType::KernelBaseType KernelType;
    KernelFactoryType *kernelFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernel = kernelFactory->CreateKernelObject(m_KernelType);
    kernel->SetKernelWidth(m_KernelWidth);
    kernel->SetSources(shapeControlPoints);
    kernel->SetWeights(momenta);
    displacementField = kernel->Convolve(shapePoints);

  } else {
    std::cerr << "Error in DisplacementFieldNormalDistribution::Sample : the dimension should be either 2 or 3."
              << std::endl;
  }

  return Superclass::m_Mean + displacementField.vectorize();
}

//template<class ScalarType>
//ScalarType
//DisplacementFieldNormalDistribution<ScalarType>
//::ComputeLogLikelihood(LinearVariableType const &obs) const {
//  const unsigned int N = Superclass::m_Mean.size();
//  assert(N == obs.n_elem());
//
//  VectorType aux = obs.vectorize();
//
//  return -0.5 * m_VarianceInverse * (Superclass::m_Mean - aux).sum_of_squares() - 0.5 * N * m_VarianceLog;
//}
//
//template<class ScalarType>
//ScalarType
//DisplacementFieldNormalDistribution<ScalarType>
//::ComputeLogLikelihood(VectorType const &obs) const {
//  const unsigned int N = Superclass::m_Mean.size();
//  assert(N == obs.size());
//
//  return -0.5 * m_VarianceInverse * (Superclass::m_Mean - obs).sum_of_squares() - 0.5 * N * m_VarianceLog;
//}


template class DisplacementFieldNormalDistribution<double>;
