/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include <itkDerivativeImageFilter.h>
#include "Diffeos.h"

/// For bug-tracking.
#include "MatrixDLM.h"

using namespace def::algebra;


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
Diffeos<ScalarType, Dimension>
::Diffeos() : Superclass(), m_T0(0.0), m_TN(1.0), m_NumberOfTimePoints(10), m_KernelType(null),
              m_KernelWidth(1.0), m_UseImprovedEuler(true), m_PaddingFactor(0.0), m_OutOfBox(true),
              m_ComputeTrueInverseFlow(false), m_UseImplicitEuler(false), m_RegressionFlag(false), m_UseFastConvolutions(false) {
  this->SetDiffeosType();
}

template<class ScalarType, unsigned int Dimension>
Diffeos<ScalarType, Dimension>
::Diffeos(const Diffeos &other) : Superclass(other) {
  this->SetDiffeosType();
  m_T0 = other.m_T0;
  m_TN = other.m_TN;
  m_NumberOfTimePoints = other.m_NumberOfTimePoints;

  m_StartPositions = other.m_StartPositions;
  m_StartMomentas = other.m_StartMomentas;

  m_PositionsT = other.m_PositionsT;
  m_MomentasT = other.m_MomentasT;

  m_KernelType = other.m_KernelType;
  m_KernelWidth = other.m_KernelWidth;

  m_ComputeTrueInverseFlow = other.m_ComputeTrueInverseFlow;
  m_UseImplicitEuler = other.m_UseImplicitEuler;
  m_RegressionFlag = other.m_RegressionFlag;
  m_UseImprovedEuler = other.m_UseImprovedEuler;

  m_DataDomain = other.m_DataDomain;
  m_BoundingBox = other.m_BoundingBox;
  m_PaddingFactor = other.m_PaddingFactor;
  m_OutOfBox = other.m_OutOfBox;

  m_MapsT = other.m_MapsT;
  m_InverseMapsT = other.m_InverseMapsT;

  m_LandmarkPointsVelocity = other.m_LandmarkPointsVelocity;
  m_LandmarkPointsT = other.m_LandmarkPointsT;
  m_AdjointPosAt0 = other.m_AdjointPosAt0;
  m_AdjointMomAt0 = other.m_AdjointMomAt0;
  m_AdjointLandmarkPointsAt0 = other.m_AdjointLandmarkPointsAt0;

  m_UseFastConvolutions = other.m_UseFastConvolutions;
}

template<class ScalarType, unsigned int Dimension>
Diffeos<ScalarType, Dimension>
::~Diffeos() {}


// template <class ScalarType, unsigned int Dimension>
// void
// Diffeos<ScalarType, Dimension>
// ::CopyInformation(const Diffeos& other)
// {
// 	
// 	Superclass::CopyInformation(other);
// 	m_KernelType = other.m_KernelType;
// 	m_KernelWidth = other.m_KernelWidth;
// 
// 	m_T0 = other.m_T0;
// 	m_TN = other.m_TN;
// 	m_NumberOfTimePoints = other.m_NumberOfTimePoints;
// 	
// 
// 	m_DataDomain = other.m_DataDomain;
// 	m_PaddingFactor = other.m_PaddingFactor;
// 	m_OutOfBox = other.m_OutOfBox;
// 	
// 	m_ComputeTrueInverseFlow = other.m_ComputeTrueInverseFlow;
// 	m_UseImprovedEuler = other.m_UseImprovedEuler;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

// template<class ScalarType, unsigned int Dimension>
// void
// Diffeos<ScalarType, Dimension>
// ::ReverseFlow()
//  {
// 	int midPoint = static_cast<int>( 0.5f*(this->GetNumberOfTimePoints() + 1) );
// 
// 	for (int t = 0; t < midPoint; t++)
// 	{
// 		MatrixType auxP = m_PositionsT[t];
// 		MatrixType auxM = - m_MomentasT[t];
// 
// 		m_PositionsT[t] = m_PositionsT[this->GetNumberOfTimePoints() - t - 1];
// 		m_MomentasT[t] = - m_MomentasT[this->GetNumberOfTimePoints() - t- 1];
// 
// 		m_PositionsT[this->GetNumberOfTimePoints() - t- 1] = auxP;
// 		m_MomentasT[this->GetNumberOfTimePoints() - t - 1] = auxM;
// 	}
// 
//  }


template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::Update() {
  if (m_NumberOfTimePoints > 1) {
    if (Superclass::IsModified()) {
      InitBoundingBox();
      Shoot();
      Superclass::UnsetModified();
      Superclass::m_DeformableObjectModified = true;
    }

    if (Superclass::m_DeformableObjectModified) {
      if (Superclass::m_IsLandmarkPoints) { FlowLandmarkPointsTrajectory(); }
      if (Superclass::m_IsImagePoints) { FlowImagePointsTrajectory(); }
      Superclass::m_DeformableObjectModified = false;
    }

  } else {

    if (Superclass::IsModified()) {
      InitBoundingBox();
      m_PositionsT.resize(1);
      m_PositionsT[0] = m_StartPositions;
      m_MomentasT.resize(1);
      m_MomentasT[0] = m_StartMomentas;
      Superclass::UnsetModified();
      Superclass::m_DeformableObjectModified = true;
    }

    if (Superclass::m_DeformableObjectModified) {
      if (Superclass::m_IsLandmarkPoints) {
        m_LandmarkPointsT.resize(1);
        m_LandmarkPointsT[0] = Superclass::m_LandmarkPoints;
      }
      if (Superclass::m_IsImagePoints) {
        m_MapsT.resize(1);
        m_InverseMapsT.resize(1);
        m_MapsT[0] = Superclass::m_ImagePoints;
        m_InverseMapsT[0] = Superclass::m_ImagePoints;
      }
      Superclass::m_DeformableObjectModified = false;
    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::InitBoundingBox() {
  ScalarType offset = m_PaddingFactor * m_KernelWidth / 2;

  MatrixType boundingBox = m_DataDomain;
  boundingBox.set_column(0, boundingBox.get_column(0) - offset);
  boundingBox.set_column(1, boundingBox.get_column(1) + offset);

  m_BoundingBox = boundingBox;
  m_OutOfBox = false;
}

template<class ScalarType, unsigned int Dimension>
bool
Diffeos<ScalarType, Dimension>
::CheckBoundingBox(MatrixListType &X, int t) {
  MatrixType Xt = X[t];
  int outOfBox = 0;

  for (int d = 0; d < Dimension; d++) {
    outOfBox += (Xt.get_column(d).min_value() < m_BoundingBox(d, 0));
    outOfBox += (Xt.get_column(d).max_value() > m_BoundingBox(d, 1));
  }

  return (outOfBox > 0);
}

template<class ScalarType, unsigned int Dimension>
std::shared_ptr<typename Diffeos<ScalarType, Dimension>::DeformableMultiObjectType>
Diffeos<ScalarType, Dimension>
::GetDeformedObjectAt(unsigned int t) const {
  if (t < 0 || t > this->GetNumberOfTimePoints() - 1)
    throw std::runtime_error("In Diffeos::GetDeformedObjectAt - Time t out of range");

  if (this->IsModified() || Superclass::m_DeformableObjectModified) {
    throw std::runtime_error(
        "In Diffeos::GetDeformedObjectAt() - The Diffeos was not updated or the objects not deformed");
  }

  MatrixType LP, IP;
  if (Superclass::m_DeformableMultiObject->GetNumberOfLandmarkKindObjects())
    LP = m_LandmarkPointsT[t];

  if (Superclass::m_DeformableMultiObject->GetNumberOfImageKindObjects())
    IP = ((m_ComputeTrueInverseFlow || m_RegressionFlag) ? m_InverseMapsT[t] : m_MapsT[m_NumberOfTimePoints - t - 1]);

  std::shared_ptr<DeformableMultiObjectType> deformedObjects
      = Superclass::m_DeformableMultiObject->DeformedMultiObject(LP, IP);

  return deformedObjects;
}

template<class ScalarType, unsigned int Dimension>
MatrixType
Diffeos<ScalarType, Dimension>
::GetDeformedControlPointsAt(unsigned int t) const {
  if (t < 0 || t > this->GetNumberOfTimePoints() - 1)
    throw std::runtime_error("In Diffeos::GetDeformedControlPointsAt() - Time t out of range");

  if (this->IsModified() || Superclass::m_DeformableObjectModified) {
    throw std::runtime_error(
        "In Diffeos::GetDeformedControlPointsAt() - The Diffeos was not updated or the objects not deformed");
  }

  return m_PositionsT[t];
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::WriteFlow(const std::vector<std::string> &name, const std::vector<std::string> &extension) {
  for (unsigned int t = 0; t < m_NumberOfTimePoints; t++) {
    std::shared_ptr<DeformableMultiObjectType> deformedObjects = GetDeformedObjectAt(t);

    std::vector<std::string> fn(name.size());
    for (int i = 0; i < name.size(); i++) {
      std::ostringstream oss;
      oss << name[i] << "__t_" << t << extension[i] << std::ends;
      fn[i] = oss.str();
    }

    deformedObjects->WriteMultiObject(fn, m_LandmarkPointsVelocity);
  }
}

template<class ScalarType, unsigned int Dimension>
MatrixType
Diffeos<ScalarType, Dimension>
::SplatResidualImage(const std::shared_ptr<DeformableMultiObjectType> target) {
  this->Update();
  return Superclass::m_DeformableMultiObject
      ->SplatDifferenceImage(
          target,
          (m_ComputeTrueInverseFlow || m_RegressionFlag) ? m_InverseMapsT[m_NumberOfTimePoints - 1]
                                                         : m_MapsT[0]);
}

template<class ScalarType, unsigned int Dimension>
MatrixType
Diffeos<ScalarType, Dimension>
::SplatResidualImages(const std::vector<std::shared_ptr<DeformableMultiObjectType>> targets,
                      const std::vector<unsigned int> &timeIndices) {
  this->Update();

  MatrixType out = Superclass::m_DeformableMultiObject
      ->SplatDifferenceImage(
          targets[0],
          (m_ComputeTrueInverseFlow || m_RegressionFlag) ? m_InverseMapsT[m_NumberOfTimePoints - 1 - timeIndices[0]]
                                                         : m_MapsT[timeIndices[0]]);

  for (unsigned int t = 1; t < timeIndices.size(); ++t)
    out += Superclass::m_DeformableMultiObject
        ->SplatDifferenceImage(
            targets[t],
            (m_ComputeTrueInverseFlow || m_RegressionFlag) ? m_InverseMapsT[m_NumberOfTimePoints - 1 - timeIndices[t]]
                                                           : m_MapsT[timeIndices[t]]);

  return out;
}

template<class ScalarType, unsigned int Dimension>
MatrixListType
Diffeos<ScalarType, Dimension>
::ParallelTransport(MatrixType const &initialMomenta,
                    MatrixType const &initialControlPoints,
                    ScalarType const &initialTime,
                    std::vector<ScalarType> const &targetTimes,
                    MatrixListType &velocities) {
  /*
   * t0 is the starting point index, initialVelocity is the velocity of the vector to be transported from this time point
   * t1 is the final point, and we fill velocities with the values of the transported velocities at the different time points.
   * We also return the values of the parallel-transported momentas, attached to the control points of the diffeo.
   */

  /// Verbose option.
  bool verbose = false;

  /// Special cases, where the transport is simply the identity :
  ///    1) Nearly zero momentas yield no motion.
  ///    2) Nearly zero initial tangent vector.
  ///    3) No target time.
  ///    4) Weird number of time points.
  const unsigned int nbTargetTimes = targetTimes.size();
  if (m_StartMomentas.sum_of_squares() < 1e-20 ||
      initialMomenta.sum_of_squares() < 1e-20 ||
      nbTargetTimes == 0 ||
      m_NumberOfTimePoints <= 1) {
    if (verbose) { std::cout << "No motion detected when computing the parallel transport." << std::endl; }

    MatrixListType out(nbTargetTimes);
    for (unsigned int t = 0; t < nbTargetTimes; ++t) {
      out[t] = initialMomenta;
    }
    return out;
  }

  /// Initial time-related operations.
  const ScalarType h = (m_TN - m_T0) / (m_NumberOfTimePoints - 1); // this is dt
  const ScalarType epsilon = h;

  // Sanity checks. Could be avoided in release mode.
  assert(nbTargetTimes > 0);
  assert(m_T0 <= initialTime && initialTime <= targetTimes[0]);
  for (unsigned int k = 1; k < nbTargetTimes; ++k) { assert(targetTimes[k - 1] < targetTimes[k]); }
  assert(targetTimes[nbTargetTimes - 1] <= m_TN + 1e-10);

  // Translates the scalar initial and target times into indexes.
  ScalarType initialIndex_continuous = (initialTime - m_T0) / h;
  unsigned int initialIndex = initialIndex_continuous;
  if ((initialIndex_continuous - initialIndex) >= 0.5) { ++initialIndex; }

  std::vector<unsigned int> targetIndices(nbTargetTimes);
  for (unsigned int k = 0; k < nbTargetTimes; ++k) {
    ScalarType targetIndex_continuous = (targetTimes[k] - m_T0) / h;
    targetIndices[k] = targetIndex_continuous;
    if ((targetIndex_continuous - targetIndices[k]) >= 0.5) { ++targetIndices[k]; }
  }

  // Further initializations.
  const unsigned int finalIndex = targetIndices[nbTargetTimes - 1];
  const unsigned int numberOfSteps = finalIndex - initialIndex + 1;

  /// Last special case, where the transport is simply the identity :
  ///    5) Nearly zero transport length.
  if (finalIndex == initialIndex) {
    if (verbose) { std::cout << "No motion detected when computing the parallel transport." << std::endl; }

    MatrixListType out(nbTargetTimes);
    for (unsigned int t = 0; t < nbTargetTimes; ++t) {
      out[t] = initialMomenta;
    }
    return out;
  }

  /// Some optional printing.
  if (verbose) {
    std::cout << "Number of control points in the matching : " << initialControlPoints.rows() << "." << std::endl;
    std::cout << "Number of control points in the regression : " << m_StartPositions.rows() << "." << std::endl;
    std::cout << "epsilon : " << epsilon << std::endl;
  }

  /// Miscellaneous initializations.
  KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(GetKernelType());
  kernelObj->SetKernelWidth(GetKernelWidth());
  MatrixType initialVelocity;
  kernelObj->SetSources(initialControlPoints);
  kernelObj->SetWeights(initialMomenta);
  initialVelocity = kernelObj->Convolve(m_PositionsT[0]);
  unsigned int numCP = this->m_StartPositions.rows();

  this->InitBoundingBox();
  MatrixType convolvKInv;
  MatrixType dPos;
  MatrixListType kGradMom;
  MatrixType CP_epsi1(numCP, Dimension, 0.); // first part of RK2
  MatrixType CP_epsi2(numCP, Dimension, 0.); // second part of RK2
  MatrixType MOM_epsi1; // first part of RK2 for the momenta, needed to compute CP_epsi1
  VectorType mom_eps_pos_temp;
  VectorType mom_eps_neg_temp;
  MatrixType mom_eps_pos(numCP, Dimension, 0.);
  MatrixType mom_eps_neg(numCP, Dimension, 0.);
  MatrixType dmom_eps_pos(numCP, Dimension, 0.);
  MatrixType dmom_eps_neg(numCP, Dimension, 0.);
  MatrixType dPos2;

  MatrixListType velocities_allSteps(numberOfSteps);
  for (unsigned int t = 0; t < numberOfSteps; ++t) { velocities_allSteps[t] = initialVelocity; }
  MatrixListType parallelTransport(numberOfSteps);

  /// These quantities should be conserved during the transport.
  ScalarType initialScalarProductWV(0.);
  ScalarType initialSquaredNormW(0.);
  ScalarType scalarProductWV(0.);
  ScalarType squaredNormW(0.);
  ScalarType squaredNormV(0.);
  MatrixType toComputeConservedQuantities;
  ScalarType alpha(1.), beta(0.);

  /// Here we compute the initial values of the conserved quantities for comparisons.
  kernelObj->SetSources(m_PositionsT[initialIndex]);
  kernelObj->SetWeights(m_MomentasT[initialIndex]);
  toComputeConservedQuantities = kernelObj->Convolve(initialControlPoints);
  for (unsigned int i = 0; i < initialMomenta.rows(); ++i)
    scalarProductWV += dot_product(initialMomenta.get_row(i), toComputeConservedQuantities.get_row(i));
  if (verbose) { std::cout << "Initial scalar product with velocities :" << scalarProductWV << std::endl; }

  kernelObj->SetSources(initialControlPoints);
  kernelObj->SetWeights(initialMomenta);
  toComputeConservedQuantities = kernelObj->Convolve(initialControlPoints);
  for (unsigned int i = 0; i < initialMomenta.rows(); ++i)
    squaredNormW += dot_product(initialMomenta.get_row(i), toComputeConservedQuantities.get_row(i));
  if (verbose) { std::cout << "Initial squared norm of w :" << squaredNormW << std::endl; }

  /// Main loop
  unsigned int index = 0;
  for (unsigned int t = initialIndex; t < finalIndex; ++t, ++index) {
    if (verbose) { std::cout << "Time step : " << t << std::endl; }
    kernelObj->SetSources(m_PositionsT[t]);
    kernelObj->SetWeights(m_MomentasT[t]);
    MatrixListType CP_epsiList;
    CP_epsiList.resize(2);
    dPos = kernelObj->Convolve(m_PositionsT[t]);

    /// This will be used to get the momenta best describing the velocities field on the control points of the diffeo.
    convolvKInv = kernelObj->ConvolveInverse(velocities_allSteps[index]);

    /// If it's the first iteration, we compute the initial scalar products and norm, to later ensure conservations.
    if (t == initialIndex) {
      kernelObj->SetSources(m_PositionsT[t]);
      kernelObj->SetWeights(m_MomentasT[t]);
      toComputeConservedQuantities = kernelObj->Convolve(m_PositionsT[t]);
      for (unsigned int i = 0; i < numCP; ++i) {
        initialScalarProductWV += dot_product(convolvKInv.get_row(i), toComputeConservedQuantities.get_row(i));
      }

      if (verbose) {
        std::cout << "Initial scalar product of proposal with velocities after projection : "
                  << initialScalarProductWV << std::endl;
      }

      kernelObj->SetWeights(convolvKInv);
      toComputeConservedQuantities = kernelObj->Convolve(m_PositionsT[t]);
      for (unsigned int i = 0; i < numCP; ++i) {
        initialSquaredNormW += dot_product(convolvKInv.get_row(i), toComputeConservedQuantities.get_row(i));
      }

      if (verbose) { std::cout << "Inital squared norm of w after projection : " << initialSquaredNormW << std::endl; }
    }

    /// We check the two conservations before updating !
    if (t > initialIndex) {
      ScalarType proposalNormSquared(0.), proposalScalarProductVelocity(0.);
      kernelObj->SetSources(m_PositionsT[t]);
      kernelObj->SetWeights(convolvKInv);
      toComputeConservedQuantities = kernelObj->Convolve(m_PositionsT[t]);
      for (unsigned int i = 0; i < numCP; ++i) {
        proposalNormSquared += dot_product(convolvKInv.get_row(i), toComputeConservedQuantities.get_row(i));
      }
      if (verbose) { std::cout << "Squared norm of proposal w : " << proposalNormSquared << std::endl; }
      for (unsigned int i = 0; i < numCP; ++i) {
        proposalScalarProductVelocity +=
            dot_product(m_MomentasT.at(t).get_row(i), toComputeConservedQuantities.get_row(i));
      }

      if (verbose) {
        std::cout << "Scalar product of proposal with velocities : " << proposalScalarProductVelocity << std::endl;
      }

      alpha = std::sqrt((initialSquaredNormW * squaredNormV - initialScalarProductWV * scalarProductWV)
                            / (proposalNormSquared * squaredNormV
                                - proposalScalarProductVelocity * proposalScalarProductVelocity));
      beta = (initialScalarProductWV - alpha * proposalScalarProductVelocity) / squaredNormV;
    }
    if (std::abs(alpha - 1.) > 0.1) {
      std::cout << ">> Warning : large alpha required to enforce the conservations. "
          "Consider decreasing the size of the steps in the scheme. (alpha = " << alpha << ")" << std::endl;
    }

    if (verbose) { std::cout << "Enforcing conservation with alpha :" << alpha << " and beta :" << beta << std::endl; }
    parallelTransport[index] = alpha * convolvKInv + beta * m_MomentasT.at(t);

    /// Here we can update our values for the conserved quantities :
    kernelObj->SetWeights(parallelTransport[index]);
    toComputeConservedQuantities = kernelObj->Convolve(m_PositionsT[t]);
    scalarProductWV = 0;
    squaredNormW = 0;
    squaredNormV = 0;
    for (unsigned int i = 0; i < numCP; ++i) {
      scalarProductWV += dot_product(m_MomentasT.at(t).get_row(i), toComputeConservedQuantities.get_row(i));
    }
    if (verbose) { std::cout << "Scalar product with velocities : " << scalarProductWV << std::endl; }
    for (unsigned int i = 0; i < numCP; ++i) {
      squaredNormW += dot_product(convolvKInv.get_row(i), toComputeConservedQuantities.get_row(i));
    }
    if (verbose) { std::cout << "Squared norm of w : " << squaredNormW << std::endl; }

    kernelObj->SetWeights(m_MomentasT[t]);
    kernelObj->SetSources(m_PositionsT[t]);
    toComputeConservedQuantities = kernelObj->Convolve(m_PositionsT[t]);
    for (unsigned int i = 0; i < numCP; ++i) {
      squaredNormV += dot_product(m_MomentasT.at(t).get_row(i), toComputeConservedQuantities.get_row(i));
    }
//        std::cout << "Squared norm of v : " << squaredNormV << std::endl;

    /// Computation of the pertubated momentum vector:
    for (unsigned int i = 0; i < numCP; i++) {
      mom_eps_pos_temp = m_MomentasT[t].get_row(i) + epsilon * convolvKInv.get_row(i);
      mom_eps_pos.set_row(i, mom_eps_pos_temp);

      mom_eps_neg_temp = m_MomentasT[t].get_row(i) - epsilon * convolvKInv.get_row(i);
      mom_eps_neg.set_row(i, (mom_eps_neg_temp));
    }

    /// We compute the hamiltonian equations with these perturbed momenta.
    kernelObj->SetWeights(mom_eps_pos);
    kGradMom = kernelObj->ConvolveGradient(m_PositionsT[t]);
    for (unsigned int i = 0; i < numCP; i++) {
      dmom_eps_pos.set_row(i, kGradMom[i].transpose() * mom_eps_pos.get_row(i));
    }

    kernelObj->SetWeights(mom_eps_neg);
    kGradMom = kernelObj->ConvolveGradient(m_PositionsT[t]);
    for (unsigned int i = 0; i < numCP; i++) {
      dmom_eps_neg.set_row(i, kGradMom[i].transpose() * mom_eps_neg.get_row(i));
    }

    /// Runge Kutta 2: computation of the middle point, computed for + epsilon.
    CP_epsi1 = m_PositionsT[t] + h / 2 * (dPos + epsilon * velocities_allSteps[index]);
    MOM_epsi1 = mom_eps_pos - h / 2 * (dmom_eps_pos); //TODO : check this !
    kernelObj->SetSources(CP_epsi1);
    kernelObj->SetWeights(MOM_epsi1);
    dPos2 = kernelObj->Convolve(CP_epsi1);
    CP_epsi2 = m_PositionsT[t] + h * dPos2;
    double div = 1 / (2 * epsilon);
    CP_epsiList[0] = CP_epsi2 * div;

    /// Runge Kutta 2: computation of the middle point, computed for - epsilon.
    CP_epsi1 = m_PositionsT[t] + h / 2 * (dPos - epsilon * velocities_allSteps[index]);
    MOM_epsi1 = mom_eps_neg - h / 2 * (dmom_eps_neg);
    kernelObj->SetSources(CP_epsi1);
    kernelObj->SetWeights(MOM_epsi1);
    dPos2 = kernelObj->Convolve(CP_epsi1);
    CP_epsi2 = m_PositionsT[t] + h * dPos2;
    CP_epsiList[1] = CP_epsi2 * div;

    /// Update the transport accordingly, in the tangent space.
    velocities_allSteps[index + 1] = (CP_epsiList[0] - CP_epsiList[1]) / epsilon;
  }

  /// Last iteration (without conservation : TODO ?).
  kernelObj->SetSources(m_PositionsT[finalIndex - 1]);
  kernelObj->SetWeights(m_MomentasT[finalIndex - 1]);
  assert(index == numberOfSteps - 1);
  convolvKInv = kernelObj->ConvolveInverse(velocities_allSteps[index]);
  parallelTransport[index] = convolvKInv;

  /// Final outputs, only at required times.
  velocities.resize(nbTargetTimes);
  MatrixListType out(nbTargetTimes);
  unsigned int outIndex = 0;
  for (unsigned int t = initialIndex; t < finalIndex + 1; ++t) {
    if (t == targetIndices[outIndex]) {
      velocities[outIndex] = velocities_allSteps[t];
      out[outIndex] = parallelTransport[t];
      ++outIndex;
    }
  }
  return out;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::Shoot() {
  unsigned int numCP =
      m_StartPositions.rows(); // WHY not creating a variable numberCP, they do not change and we could set it as const

  MatrixListType &outPos = m_PositionsT;
  MatrixListType &outMoms = m_MomentasT;

  outPos.resize(m_NumberOfTimePoints);
  outMoms.resize(m_NumberOfTimePoints);
  for (unsigned int t = 0; t < m_NumberOfTimePoints; t++) {
    outPos[t] = m_StartPositions;
    outMoms[t].set_size(numCP, Dimension);
    outMoms[t].fill(0);
  }
  outPos[0] = m_StartPositions;
  outMoms[0] = m_StartMomentas;

  // Special case: nearly zero momentas yield no motion
  if (outMoms[0].frobenius_norm() < 1e-20)
    return;

  //ScalarType timeStep = 1.0 / (m_NumberOfTimePoints-1); // Why not T0 and Tm
  ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

  KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
  kernelObj->SetKernelWidth(this->GetKernelWidth());

  for (unsigned int t = 0; t < (m_NumberOfTimePoints - 1); t++) {
    kernelObj->SetSources(outPos[t]);
    kernelObj->SetWeights(outMoms[t]);

    std::vector<MatrixType> kGradMom = kernelObj->ConvolveGradient(outPos[t]);

    MatrixType dPos = kernelObj->Convolve(outPos[t]);
    MatrixType dMom(numCP, Dimension, 0);
    for (unsigned int i = 0; i < numCP; i++)
      dMom.set_row(i, (kGradMom[i].transpose() * outMoms[t].get_row(i)));

    outPos[t + 1] = outPos[t] + dPos * dt;
    outMoms[t + 1] = outMoms[t] - dMom * dt;

    //		// Heun's method
    //		if (m_UseImprovedEuler)
    //		{
    //			kernelObj->SetSources(outPos[t+1]);
    //			kernelObj->SetWeights(outMoms[t+1]);
    //
    //			kGradMom = kernelObj->ConvolveGradient(outPos[t+1]);
    //
    //			typename KernelType::MatrixType dPos2 =
    //					kernelObj->Convolve(outPos[t+1]);
    //
    //			typename KernelType::MatrixType dMom2(numCP, Dimension, 0);
    //			for (unsigned int i = 0; i < numCP; i++)
    //				dMom2.set_row(i,
    //						(kGradMom[i].transpose() * outMoms[t+1].get_row(i)) );
    //
    //			outPos[t+1] = outPos[t] + (dPos + dPos2) * (timeStep * 0.5f);
    //			outMoms[t+1] = outMoms[t] - (dMom + dMom2) * (timeStep * 0.5f);
    //		}

//    if (this->CheckBoundingBox(outPos, t + 1)) {
//      std::cout << ">> Geodesic shooting : out of box at time t = " << t + 1 << "." << std::endl;
//      break;
//    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::FlowLandmarkPointsTrajectory() {

  //ScalarType dt = 1.0 / (m_NumberOfTimePoints - 1);
  ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  typedef typename KernelFactoryType::KernelBaseType KernelType;
  KernelFactoryType *kFactory = KernelFactoryType::Instantiate();

  std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
  kernelObj->SetKernelWidth(this->GetKernelWidth());

  MatrixListType posT = this->GetTrajectoryPositions();
  MatrixListType momT = this->GetTrajectoryMomentas();

  m_LandmarkPointsT.resize(m_NumberOfTimePoints);
  m_LandmarkPointsVelocity.resize(m_NumberOfTimePoints);
  for (unsigned int t = 0; t < m_NumberOfTimePoints; t++) {
    m_LandmarkPointsT[t] = Superclass::m_LandmarkPoints;
    m_LandmarkPointsVelocity[t].fill(0.0);
  }

  if (momT[0].frobenius_norm() < 1e-20)
    return;

  for (unsigned int t = 0; t < m_NumberOfTimePoints - 1; t++) {
    kernelObj->SetSources(posT[t]);
    kernelObj->SetWeights(momT[t]);

    m_LandmarkPointsVelocity[t] = kernelObj->Convolve(m_LandmarkPointsT[t]);
    m_LandmarkPointsT[t + 1] = m_LandmarkPointsT[t] + (m_LandmarkPointsVelocity[t] * dt);

    if (this->ImprovedEuler()) {
      kernelObj->SetSources(posT[t + 1]);
      kernelObj->SetWeights(momT[t + 1]);

      m_LandmarkPointsVelocity[t + 1] = kernelObj->Convolve(m_LandmarkPointsT[t + 1]);
      m_LandmarkPointsT[t + 1] =
          m_LandmarkPointsT[t] + (m_LandmarkPointsVelocity[t] + m_LandmarkPointsVelocity[t + 1]) * (dt * 0.5f);
    }

//    if (this->CheckBoundingBox(m_LandmarkPointsT, t + 1)) {
//      std::cout << ">> Landmark deformation: out of box at time t = " << t + 1 << "." << std::endl;
//      break;
//    }
  }
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::FlowImagePointsTrajectory() {

  if (m_ComputeTrueInverseFlow) { IntegrateImagePointsWithTrueInverseFlow(); }
  else { IntegrateImagePointsBackward(); }
}

// Computes \phi_t\circ\phi_1^{-1}, which is equal to \phi_1^{-1} for t = 0
template<class ScalarType, unsigned int Dimension>
void Diffeos<ScalarType, Dimension>
::IntegrateImagePointsBackward() {
  if (!m_RegressionFlag) {
    m_MapsT.resize(m_NumberOfTimePoints);
    for (unsigned int t = 0; t < m_NumberOfTimePoints; ++t) { m_MapsT[t] = Superclass::m_ImagePoints; }

    /// Special case: nearly zero momentas yield no motion
    if (m_MomentasT[0].frobenius_norm() < 1e-20) { return; }

    /// Initializes the Euler time step.
    ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

    /// Initializes the kernel object.
    KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
    kernelObj->SetKernelWidth(m_KernelWidth);

    /// Backward integration, Euler scheme.
    for (unsigned int t = this->GetNumberOfTimePoints() - 1; t >= 1; --t) {
      kernelObj->SetSources(m_PositionsT[t]);
      kernelObj->SetWeights(m_MomentasT[t]);

      MatrixType dY;
      if (m_UseFastConvolutions) {
        dY = kernelObj->ConvolveImageFast(m_MapsT[t], Superclass::m_DownSampledImage);
      }
      else { dY = kernelObj->Convolve(m_MapsT[t]); }

      m_MapsT[t - 1] = m_MapsT[t] - dY * dt;

      /// Heun's method.
      if (m_UseImprovedEuler) {
        kernelObj->SetSources(m_PositionsT[t - 1]);
        kernelObj->SetWeights(m_MomentasT[t - 1]);

        MatrixType dY2;
        if (m_UseFastConvolutions) {
          dY2 = kernelObj->ConvolveImageFast(m_MapsT[t - 1], Superclass::m_DownSampledImage);
        }
        else { dY2 = kernelObj->Convolve(m_MapsT[t - 1]); }

        m_MapsT[t - 1] = m_MapsT[t] - (dY + dY2) * (dt * 0.5f);
      }

      if (this->CheckBoundingBox(m_MapsT, t - 1)) {
        std::cout << ">> Image deformation: out of box at time t = " << t - 1 << "." << std::endl;
        return;
      }
    }

  } else {
    /// Initialization of the inverse maps.
    m_InverseMapsT.resize(m_NumberOfTimePoints);
    for (unsigned int t = 0; t < m_NumberOfTimePoints; ++t) { m_InverseMapsT[t] = Superclass::m_ImagePoints; }

    /// Special case: nearly zero momentas yield no motion
    if (m_MomentasT[0].frobenius_norm() < 1e-20) { return; }

    /// Initializes the Euler time step.
    ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

    /// Initializes the kernel object.
    KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
    std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
    kernelObj->SetKernelWidth(m_KernelWidth);

    for (unsigned int tt = 1; tt < m_NumberOfTimePoints; ++tt) {
      /// Backward integration, Euler scheme.
      MatrixType aux = Superclass::m_ImagePoints;
      for (unsigned int t = tt; t > 0; --t) {
        kernelObj->SetSources(m_PositionsT[t]);
        kernelObj->SetWeights(m_MomentasT[t]);

        MatrixType dY;
        if (m_UseFastConvolutions) {
          dY = kernelObj->ConvolveImageFast(aux, Superclass::m_DownSampledImage);
        }
        else { dY = kernelObj->Convolve(aux); }

        /// Heun's method.
        if (m_UseImprovedEuler) {
          kernelObj->SetSources(m_PositionsT[t - 1]);
          kernelObj->SetWeights(m_MomentasT[t - 1]);

          MatrixType dY2;
          if (m_UseFastConvolutions) {
            dY2 = kernelObj->ConvolveImageFast(aux - dY * dt, Superclass::m_DownSampledImage);
          }
          else { dY2 = kernelObj->Convolve(aux - dY * dt); }

          aux = aux - (dY + dY2) * (dt * 0.5f);
        } else {
          aux = aux - dY * dt;
        }
      }
      m_InverseMapsT[tt] = aux;
    }
  }
}


template<class ScalarType, unsigned int Dimension>
MatrixType
Diffeos<ScalarType, Dimension>
::ComputeExplicitEulerStep(const unsigned int t,
                           const MatrixType &v,
                           const std::vector<unsigned int> &size) {
  unsigned int spatialApproximationOrder = 1;

  switch (spatialApproximationOrder) {
    case 1 : return ComputeExplicitEulerStep_1(t, v, size);
    case 4 : return ComputeExplicitEulerStep_4(t, v, size);
    default : return ComputeExplicitEulerStep_1(t, v, size);
  }
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::SolveImplicitEulerStep(const unsigned int t,
                         const MatrixType &v,
                         const std::vector<unsigned int> &size) {
  unsigned int spatialApproximationOrder = 2;

  switch (spatialApproximationOrder) {
    case 2 : return SolveImplicitEulerStep_2(t, v, size);
    default : return SolveImplicitEulerStep_2(t, v, size);
  }
}

// At order 1 in space.
template<class ScalarType, unsigned int Dimension>
MatrixType
Diffeos<ScalarType, Dimension>
::ComputeExplicitEulerStep_1(const unsigned int t,
                             const MatrixType &v,
                             const std::vector<unsigned int> &size) {
  MatrixType dY;

  if (Dimension == 2) {
    dY.set_size(size[0] * size[1], 2);
    dY.fill(0.0);
    unsigned int p = 0;

    /// X direction.
    for (unsigned int j = 0; j < size[1]; ++j) {
      // Top, i = 0 (forward).
      p = j * size[0];
      dY.increment_row(p, -v(p, 0) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
      // Core (central).
      for (unsigned int i = 1; i < size[0] - 1; ++i) {
        p = i + j * size[0];
        dY.increment_row(p, -0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
        dY.increment_row(p, 0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
      }
      // Bottom, i = size[0] - 1 (backward).
      p = (j + 1) * size[0] - 1;
      dY.increment_row(p, v(p, 0) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, -v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
    }

    /// Y direction.
    for (unsigned int i = 0; i < size[0]; ++i) {
      // Top, j = 0 (forward).
      p = i;
      dY.increment_row(p, -v(p, 1) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
      // Core (central).
      for (unsigned int j = 1; j < size[1] - 1; ++j) {
        p = i + j * size[0];
        dY.increment_row(p, -0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
        dY.increment_row(p, 0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
      }
      // Bottom, j = size[1] - 1 (backward).
      p = i + (size[1] - 1) * size[0];
      dY.increment_row(p, v(p, 1) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, -v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
    }
  } else if (Dimension == 3) {
    dY.set_size(size[0] * size[1] * size[2], 3);
    dY.fill(0.0);

    unsigned int p = 0;
    for (unsigned int k = 0; k < size[2]; ++k) {
      /// X direction.
      for (unsigned int j = 0; j < size[1]; ++j) {
        // Top, i = 0 (forward).
        p = j * size[0] + k * size[0] * size[1];
        dY.increment_row(p, -v(p, 0) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
        // Core (central).
        for (unsigned int i = 1; i < size[0] - 1; ++i) {
          p = i + j * size[0] + k * size[0] * size[1];
          dY.increment_row(p, -0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
          dY.increment_row(p, 0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
        }
        // Bottom, i = size[0] - 1 (backward).
        p = (j + 1) * size[0] - 1 + k * size[0] * size[1];
        dY.increment_row(p, v(p, 0) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
      }

      /// Y direction.
      for (unsigned int i = 0; i < size[0]; ++i) {
        // Top, j = 0 (forward).
        p = i + k * size[0] * size[1];
        dY.increment_row(p, -v(p, 1) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
        // Core (central).
        for (unsigned int j = 1; j < size[1] - 1; ++j) {
          p = i + j * size[0] + k * size[0] * size[1];
          dY.increment_row(p, -0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
          dY.increment_row(p, 0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
        }
        // Bottom, j = size[1] - 1 (backward).
        p = i + (size[1] - 1) * size[0] + k * size[0] * size[1];
        dY.increment_row(p, v(p, 1) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
      }
    }

    /// Z direction.
    for (unsigned int i = 0; i < size[0]; ++i) {
      for (unsigned int j = 0; j < size[1]; ++j) {
        // Top, k = 0 (forward).
        p = i + j * size[0];
        dY.increment_row(p, -v(p, 2) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, v(p, 2) * m_InverseMapsT[t].get_row(p + size[0] * size[1]));
        // Core (central).
        for (unsigned int k = 1; k < size[2] - 1; ++k) {
          p = i + j * size[0] + k * size[0] * size[1];
          dY.increment_row(p, -0.5 * v(p, 2) * m_InverseMapsT[t].get_row(p - size[0] * size[1]));
          dY.increment_row(p, 0.5 * v(p, 2) * m_InverseMapsT[t].get_row(p + size[0] * size[1]));
        }
        // Bottom, k = size[2] - 1 (backward).
        p = i + j * size[0] + (size[2] - 1) * size[0] * size[1];
        dY.increment_row(p, v(p, 2) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -v(p, 2) * m_InverseMapsT[t].get_row(p - size[0] * size[1]));
      }
    }
  } else {
    std::cerr << "Error : in Diffeos::ComputeExplicitEulerStep, the dimension must be either 2 or 3."
              << std::endl;
  }

  return dY;
}

// At order 4 in space.
template<class ScalarType, unsigned int Dimension>
MatrixType
Diffeos<ScalarType, Dimension>
::ComputeExplicitEulerStep_4(const unsigned int t,
                             const MatrixType &v,
                             const std::vector<unsigned int> &size) {
  MatrixType dY;

  if (Dimension == 2) {
    dY.set_size(size[0] * size[1], 2);
    dY.fill(0.0);
    unsigned int p = 0;

    /// X direction.
    for (unsigned int j = 0; j < size[1]; ++j) {
      // Top 1, i = 0 (forward).
      p = j * size[0];
      dY.increment_row(p, -25. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, 4 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
      dY.increment_row(p, -3 * v(p, 0) * m_InverseMapsT[t].get_row(p + 2));
      dY.increment_row(p, 4. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p + 3));
      dY.increment_row(p, -0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p + 4));
      // Top 2, i = 1 (forward).
      p = 1 + j * size[0];
      dY.increment_row(p, -0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
      dY.increment_row(p, -5. / 6 * v(p, 0) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, 1.5 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
      dY.increment_row(p, -0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p + 2));
      dY.increment_row(p, 1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p + 3));
      // Core (central).
      for (unsigned int i = 2; i < size[0] - 2; ++i) {
        p = i + j * size[0];
        dY.increment_row(p, 1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p - 2));
        dY.increment_row(p, -2. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
        dY.increment_row(p, 2. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
        dY.increment_row(p, -1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p + 2));
      }
      // Bottom 1, i = size[0] - 2 (backward).
      p = (j + 1) * size[0] - 2;
      dY.increment_row(p, +0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
      dY.increment_row(p, +5. / 6 * v(p, 0) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, -1.5 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
      dY.increment_row(p, +0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p - 2));
      dY.increment_row(p, -1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p - 3));
      // Bottom 2, i = size[0] - 1 (backward).
      p = (j + 1) * size[0] - 1;
      dY.increment_row(p, 25. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, -4 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
      dY.increment_row(p, 3 * v(p, 0) * m_InverseMapsT[t].get_row(p - 2));
      dY.increment_row(p, -4. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p - 3));
      dY.increment_row(p, 0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p - 4));
    }

    /// Y direction.
    for (unsigned int i = 0; i < size[0]; ++i) {
      // Top 1, j = 0 (forward).
      p = i;
      dY.increment_row(p, -25. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, 4 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
      dY.increment_row(p, -3 * v(p, 1) * m_InverseMapsT[t].get_row(p + 2 * size[0]));
      dY.increment_row(p, 4. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p + 3 * size[0]));
      dY.increment_row(p, -0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p + 4 * size[0]));
      // Top 2, j = 1 (forward).
      p = i + size[0];
      dY.increment_row(p, -0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
      dY.increment_row(p, -5. / 6 * v(p, 1) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, 1.5 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
      dY.increment_row(p, -0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p + 2 * size[0]));
      dY.increment_row(p, 1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p + 3 * size[0]));
      // Core (central).
      for (unsigned int j = 2; j < size[1] - 2; ++j) {
        p = i + j * size[0];
        dY.increment_row(p, 1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p - 2 * size[0]));
        dY.increment_row(p, -2. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
        dY.increment_row(p, 2. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
        dY.increment_row(p, -1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p + 2 * size[0]));
      }
      // Bottom 1, j = size[1] - 2 (backward).
      p = i + (size[1] - 2) * size[0];
      dY.increment_row(p, +0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
      dY.increment_row(p, +5. / 6 * v(p, 1) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, -1.5 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
      dY.increment_row(p, +0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p - 2 * size[0]));
      dY.increment_row(p, -1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p - 3 * size[0]));
      // Bottom 2, j = size[1] - 1 (backward).
      p = i + (size[1] - 1) * size[0];
      dY.increment_row(p, 25. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p));
      dY.increment_row(p, -4 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
      dY.increment_row(p, 3 * v(p, 1) * m_InverseMapsT[t].get_row(p - 2 * size[0]));
      dY.increment_row(p, -4. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p - 3 * size[0]));
      dY.increment_row(p, 0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p - 4 * size[0]));
    }
  } else if (Dimension == 3) {
    dY.set_size(size[0] * size[1] * size[2], 3);
    dY.fill(0.0);

    unsigned int p = 0;
    for (unsigned int k = 0; k < size[2]; ++k) {
      /// X direction.
      for (unsigned int j = 0; j < size[1]; ++j) {
        // Top 1, i = 0 (forward).
        p = j * size[0] + k * size[0] * size[1];
        dY.increment_row(p, -25. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, 4 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
        dY.increment_row(p, -3 * v(p, 0) * m_InverseMapsT[t].get_row(p + 2));
        dY.increment_row(p, 4. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p + 3));
        dY.increment_row(p, -0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p + 4));
        // Top 2, i = 1 (forward).
        p = 1 + j * size[0] + k * size[0] * size[1];
        dY.increment_row(p, -0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
        dY.increment_row(p, -5. / 6 * v(p, 0) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, 1.5 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
        dY.increment_row(p, -0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p + 2));
        dY.increment_row(p, 1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p + 3));
        // Core (central).
        for (unsigned int i = 2; i < size[0] - 2; ++i) {
          p = i + j * size[0] + k * size[0] * size[1];
          dY.increment_row(p, 1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p - 2));
          dY.increment_row(p, -2. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
          dY.increment_row(p, 2. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
          dY.increment_row(p, -1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p + 2));
        }
        // Bottom 1, i = size[0] - 2 (backward).
        p = (j + 1) * size[0] - 2 + k * size[0] * size[1];
        dY.increment_row(p, +0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p + 1));
        dY.increment_row(p, +5. / 6 * v(p, 0) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -1.5 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
        dY.increment_row(p, +0.5 * v(p, 0) * m_InverseMapsT[t].get_row(p - 2));
        dY.increment_row(p, -1. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p - 3));
        // Bottom 2, i = size[0] - 1 (backward).
        p = (j + 1) * size[0] - 1 + k * size[0] * size[1];
        dY.increment_row(p, 25. / 12 * v(p, 0) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -4 * v(p, 0) * m_InverseMapsT[t].get_row(p - 1));
        dY.increment_row(p, 3 * v(p, 0) * m_InverseMapsT[t].get_row(p - 2));
        dY.increment_row(p, -4. / 3 * v(p, 0) * m_InverseMapsT[t].get_row(p - 3));
        dY.increment_row(p, 0.25 * v(p, 0) * m_InverseMapsT[t].get_row(p - 4));
      }

      /// Y direction.
      for (unsigned int i = 0; i < size[0]; ++i) {
        // Top 1, j = 0 (forward).
        p = i + k * size[0] * size[1];
        dY.increment_row(p, -25. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, 4 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
        dY.increment_row(p, -3 * v(p, 1) * m_InverseMapsT[t].get_row(p + 2 * size[0]));
        dY.increment_row(p, 4. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p + 3 * size[0]));
        dY.increment_row(p, -0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p + 4 * size[0]));
        // Top 2, j = 1 (forward).
        p = i + size[0] + k * size[0] * size[1];
        dY.increment_row(p, -0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
        dY.increment_row(p, -5. / 6 * v(p, 1) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, 1.5 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
        dY.increment_row(p, -0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p + 2 * size[0]));
        dY.increment_row(p, 1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p + 3 * size[0]));
        // Core (central).
        for (unsigned int j = 2; j < size[1] - 2; ++j) {
          p = i + j * size[0] + k * size[0] * size[1];
          dY.increment_row(p, 1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p - 2 * size[0]));
          dY.increment_row(p, -2. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
          dY.increment_row(p, 2. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
          dY.increment_row(p, -1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p + 2 * size[0]));
        }
        // Bottom 1, j = size[1] - 2 (backward).
        p = i + (size[1] - 2) * size[0] + k * size[0] * size[1];
        dY.increment_row(p, +0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p + size[0]));
        dY.increment_row(p, +5. / 6 * v(p, 1) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -1.5 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
        dY.increment_row(p, +0.5 * v(p, 1) * m_InverseMapsT[t].get_row(p - 2 * size[0]));
        dY.increment_row(p, -1. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p - 3 * size[0]));
        // Bottom 2, j = size[1] - 1 (backward).
        p = i + (size[1] - 1) * size[0] + k * size[0] * size[1];
        dY.increment_row(p, 25. / 12 * v(p, 1) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -4 * v(p, 1) * m_InverseMapsT[t].get_row(p - size[0]));
        dY.increment_row(p, 3 * v(p, 1) * m_InverseMapsT[t].get_row(p - 2 * size[0]));
        dY.increment_row(p, -4. / 3 * v(p, 1) * m_InverseMapsT[t].get_row(p - 3 * size[0]));
        dY.increment_row(p, 0.25 * v(p, 1) * m_InverseMapsT[t].get_row(p - 4 * size[0]));
      }
    }

    /// Z direction.
    for (unsigned int i = 0; i < size[0]; ++i) {
      for (unsigned int j = 0; j < size[1]; ++j) {
        // Top 1, k = 0 (forward).
        p = i + j * size[0];
        dY.increment_row(p, -25. / 12 * v(p, 2) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, 4 * v(p, 2) * m_InverseMapsT[t].get_row(p + size[0] * size[1]));
        dY.increment_row(p, -3 * v(p, 2) * m_InverseMapsT[t].get_row(p + 2 * size[0] * size[1]));
        dY.increment_row(p, 4. / 3 * v(p, 2) * m_InverseMapsT[t].get_row(p + 3 * size[0] * size[1]));
        dY.increment_row(p, -0.25 * v(p, 2) * m_InverseMapsT[t].get_row(p + 4 * size[0] * size[1]));
        // Top 2, k = 1 (forward).
        p = i + j * size[0] + size[0] * size[1];
        dY.increment_row(p, -0.25 * v(p, 2) * m_InverseMapsT[t].get_row(p - size[0] * size[1]));
        dY.increment_row(p, -5. / 6 * v(p, 2) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, 1.5 * v(p, 2) * m_InverseMapsT[t].get_row(p + size[0] * size[1]));
        dY.increment_row(p, -0.5 * v(p, 2) * m_InverseMapsT[t].get_row(p + 2 * size[0] * size[1]));
        dY.increment_row(p, 1. / 12 * v(p, 2) * m_InverseMapsT[t].get_row(p + 3 * size[0] * size[1]));
        // Core (central).
        for (unsigned int k = 2; k < size[2] - 2; ++k) {
          p = i + j * size[0] + k * size[0] * size[1];
          dY.increment_row(p, 1. / 12 * v(p, 2) * m_InverseMapsT[t].get_row(p - 2 * size[0] * size[1]));
          dY.increment_row(p, -2. / 3 * v(p, 2) * m_InverseMapsT[t].get_row(p - size[0] * size[1]));
          dY.increment_row(p, 2. / 3 * v(p, 2) * m_InverseMapsT[t].get_row(p + size[0] * size[1]));
          dY.increment_row(p, -1. / 12 * v(p, 2) * m_InverseMapsT[t].get_row(p + 2 * size[0] * size[1]));
        }
        // Bottom 1, k = size[2] - 2 (backward).
        p = i + j * size[0] + (size[2] - 2) * size[0] * size[1];
        dY.increment_row(p, +0.25 * v(p, 2) * m_InverseMapsT[t].get_row(p + size[0] * size[1]));
        dY.increment_row(p, +5. / 6 * v(p, 2) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -1.5 * v(p, 2) * m_InverseMapsT[t].get_row(p - size[0] * size[1]));
        dY.increment_row(p, +0.5 * v(p, 2) * m_InverseMapsT[t].get_row(p - 2 * size[0] * size[1]));
        dY.increment_row(p, -1. / 12 * v(p, 2) * m_InverseMapsT[t].get_row(p - 3 * size[0] * size[1]));
        // Bottom 2, k = size[2] - 1 (backward).
        p = i + j * size[0] + (size[2] - 1) * size[0] * size[1];
        dY.increment_row(p, 25. / 12 * v(p, 2) * m_InverseMapsT[t].get_row(p));
        dY.increment_row(p, -4 * v(p, 2) * m_InverseMapsT[t].get_row(p - size[0] * size[1]));
        dY.increment_row(p, 3 * v(p, 2) * m_InverseMapsT[t].get_row(p - 2 * size[0] * size[1]));
        dY.increment_row(p, -4. / 3 * v(p, 2) * m_InverseMapsT[t].get_row(p - 3 * size[0] * size[1]));
        dY.increment_row(p, 0.25 * v(p, 2) * m_InverseMapsT[t].get_row(p - 4 * size[0] * size[1]));
      }
    }
  } else {
    std::cerr << "Error : in Diffeos::ComputeExplicitEulerStep, the dimension must be either 2 or 3."
              << std::endl;
  }

  return dY;
}

// At order 2 in space.
template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::SolveImplicitEulerStep_2(const unsigned int t,
                           const MatrixType &v,
                           const std::vector<unsigned int> &size) {
  if (Dimension == 2) {
    MatrixType M(size[0] * size[1], size[0] * size[1], 0.0);
    unsigned int p = 0;

    /// X direction.
    for (unsigned int j = 0; j < size[1]; ++j) {
      // Top, i = 0 (forward).
      p = j * size[0];
      M(p, p) += -1.5 * v(p, 0);
      M(p, p + 1) += 2 * v(p, 0);
      M(p, p + 2) += -0.5 * v(p, 0);
      for (unsigned int i = 1; i < size[0] - 1; ++i) {
        p = i + j * size[0];
        M(p, p - 1) += -0.5 * v(p, 0);
        M(p, p + 1) += 0.5 * v(p, 0);
      }
      // Bottom, i = size[0] - 1 (backward).
      p = (j + 1) * size[0] - 1;
      M(p, p) += 1.5 * v(p, 0);
      M(p, p - 1) += -2 * v(p, 0);
      M(p, p - 2) += 0.5 * v(p, 0);
    }

    /// Y direction.
    for (unsigned int i = 0; i < size[0]; ++i) {
      // Top, j = 0 (forward).
      p = i;
      M(p, p) += -1.5 * v(p, 1);
      M(p, p + size[0]) += 2 * v(p, 1);
      M(p, p + 2 * size[0]) += -0.5 * v(p, 1);
      for (unsigned int j = 1; j < size[0] - 1; ++j) {
        p = i + j * size[0];
        M(p, p - size[0]) += -0.5 * v(p, 1);
        M(p, p + size[0]) += 0.5 * v(p, 1);
      }
      // Bottom, j = size[1] - 1 (backward).
      p = i + (size[1] - 1) * size[0];
      M(p, p) += 1.5 * v(p, 1);
      M(p, p - size[0]) += -2 * v(p, 1);
      M(p, p - 2 * size[0]) += 0.5 * v(p, 1);
    }

    m_InverseMapsT[t + 1] = solve(
        diagonal_matrix<ScalarType>(size[0] * size[1], 1.0) + M, m_InverseMapsT[t]);
  } else if (Dimension == 3) {
    std::cerr << "Error : in Diffeos::ComputeDiscretizedEquations, the dimension == 3 case is not available yet."
              << std::endl;
  } else
    std::cerr << "Error : in Diffeos::ComputeDiscretizedEquations, the dimension must be either 2 or 3."
              << std::endl;
}

// Computes the flow \phi_t^{-1}
template<class ScalarType, unsigned int Dimension>
void Diffeos<ScalarType, Dimension>
::IntegrateImagePointsWithTrueInverseFlow() {
  /// Initialize all the time points to be the image at the first time point.
  m_InverseMapsT.resize(m_NumberOfTimePoints);
  for (long t = 0; t < m_NumberOfTimePoints; t++)
    m_InverseMapsT[t] = Superclass::m_ImagePoints;

  /// Special case: nearly zero momentas yield no motion.
  if (m_MomentasT[0].frobenius_norm() < 1e-20)
    return;

  /// The time step.
  ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

  // A list of Dimension images, each image stores a spatial coordinate representing the physical location of the image points
  std::vector<ImagePointerType> Yt;
  Yt.resize(Dimension);
  for (unsigned int d = 0; d < Dimension; d++) {
    // Create an image from each Dimension
    ImagePointerType
        img = GridFunctionsType::VectorToImage(Superclass::m_DownSampledImage, Superclass::m_ImagePoints.get_column(d));
    // ImagePointerType img = GridFunctionsType::VectorToImage(downSampledWorkingImage, downSampledY1.get_column(d));
    Yt[d] = img;
  }

  // The kernel is for computing v_t(y0)
  KernelFactoryType *kFactory = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(this->GetKernelType());
  kernelObj->SetKernelWidth(this->GetKernelWidth());

  /// Loop over time for the Euler scheme.
  for (unsigned long t = 0; t < m_NumberOfTimePoints - 1; ++t) {
    kernelObj->SetSources(m_PositionsT[t]);
    kernelObj->SetWeights(m_MomentasT[t]);
    MatrixType VtY0;
    if (m_UseFastConvolutions)
      VtY0 = kernelObj->ConvolveImageFast(Superclass::m_ImagePoints, Superclass::m_DownSampledImage);
    else
      VtY0 = kernelObj->Convolve(Superclass::m_ImagePoints);

    // Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
    std::vector<ImagePointerType> gradImages;
    gradImages.resize(Dimension * Dimension);
    unsigned int indx = 0;
    // Compute the gradient in all directions
    for (unsigned int dim1 = 0; dim1 < Dimension; dim1++) {
      for (unsigned int dim2 = 0; dim2 < Dimension; dim2++) {
        typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
        typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
        derivf->SetInput(Yt[dim1]);
        derivf->SetDirection(dim2);
        derivf->SetOrder(1);
        derivf->SetUseImageSpacingOn();
        derivf->Update();

        gradImages[indx] = derivf->GetOutput();
        indx++;
      }
    }

    // Now we have all the information we need to compute dY, but we need to loop over all the pixels
    // and construct the 3x3 jacobian matrix and 3x1 v_t(y0)
    ImageIteratorType it(Superclass::m_DownSampledImage, Superclass::m_DownSampledImage->GetLargestPossibleRegion());
    MatrixType dY(Superclass::m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels(), Dimension);

    // Loop over the grid to construct jacobian matrices and compute dY
    unsigned int k = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
      typedef typename ImageType::IndexType ImageIndexType;
      ImageIndexType ind = it.GetIndex();
      MatrixType jacobian(Dimension, Dimension);

      unsigned int indx = 0;
      // Build the (DimensionxDimension) jacobian matrix
      for (unsigned int dim1 = 0; dim1 < Dimension; dim1++) {
        for (unsigned int dim2 = 0; dim2 < Dimension; dim2++) {
          ImagePointerType curGradImage = gradImages[indx];
          jacobian(dim1, dim2) = curGradImage->GetPixel(ind);
          indx++;
        }
      }

      // Get the (Dimension) vector corresponding to this points v_t(y0)
      VectorType VtY0k = VtY0.get_row(k);

      // dY = d_y0 \phi_t * v_t(y0)
      dY.set_row(k, jacobian * VtY0k);
      k++;
    }

    /// Updates the inverse map by Euler scheme.
    m_InverseMapsT[t + 1] = m_InverseMapsT[t] - dY * dt;

    // Update Yt using the m_InverseMap matrix we just updated using Euler
    MatrixType inversem1 = m_InverseMapsT[t + 1];
    for (unsigned int d = 0; d < Dimension; d++) {
      // Create an image from each Dimension
      ImagePointerType img = GridFunctionsType::VectorToImage(Superclass::m_DownSampledImage, inversem1.get_column(d));
      Yt[d] = img;
    }

//		if (this->CheckBoundingBox(m_InverseMapsT, t+1))
//		{
//			std::cout << "Image deformation: out of box at time t = " << t+1 << std::endl;
//			return;
//		}

  }
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::IntegrateAdjointEquations(MatrixType &InitialConditionsLandmarkPoints,
                            MatrixType &InitialConditionsImagePoints) {
  MatrixListType ListInitialConditionsLandmarkPoints;
  MatrixListType ListInitialConditionsImagePoints;
  ListInitialConditionsLandmarkPoints.resize(1);
  ListInitialConditionsImagePoints.resize(1);

  ListInitialConditionsLandmarkPoints[0] = InitialConditionsLandmarkPoints;
  ListInitialConditionsImagePoints[0] = InitialConditionsImagePoints;

  std::vector<unsigned int> times;
  times.resize(0);
  this->IntegrateAdjointEquations(ListInitialConditionsLandmarkPoints, ListInitialConditionsImagePoints, times);
}

template<class ScalarType, unsigned int Dimension>
void
Diffeos<ScalarType, Dimension>
::IntegrateAdjointEquations(MatrixListType &InitialConditionsLandmarkPoints,
                            MatrixListType &InitialConditionsImagePoints,
                            std::vector<unsigned int> jumpTimes) {
  // Upsample image maps, since initial condition of image objects are at full resolution
  // ImagePointerType image = Superclass::m_Template->GetTemplateObjects()->GetImage();
  // The path of the image points over time
  MatrixListType fullResMapsT(m_NumberOfTimePoints);
  typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;

  if (Superclass::m_IsImagePoints) {
    for (unsigned int t = 0; t < m_NumberOfTimePoints; t++)
      fullResMapsT[t] = GridFunctionsType::UpsampleImagePoints(
          Superclass::m_Image, Superclass::m_DownSampledImage,
          (m_ComputeTrueInverseFlow || m_RegressionFlag) ? m_InverseMapsT[t] : m_MapsT[t]);
  }

  typedef AdjointEquationsIntegrator<ScalarType, Dimension> AEIntegrator;
  AEIntegrator *integrator = new AEIntegrator();
  integrator->SetControlPointsTrajectory(m_PositionsT);
  integrator->SetMomentaTrajectory(m_MomentasT);
  if (Superclass::m_IsLandmarkPoints) {
    integrator->SetLandmarkPointsTrajectory(m_LandmarkPointsT);
    integrator->SetInitialConditionsLandmarkPoints(InitialConditionsLandmarkPoints);
  }
  if (Superclass::m_IsImagePoints) {
    integrator->SetImagePointsTrajectory(fullResMapsT);
    integrator->SetInitialConditionsImagePoints(InitialConditionsImagePoints);
  }
  integrator->SetKernelWidth(m_KernelWidth);
  integrator->SetKernelType(m_KernelType);
  integrator->SetNumberOfTimePoints(m_NumberOfTimePoints);
  integrator->SetFullResolutionImage(Superclass::m_Image);
  integrator->SetDownSampledImage(Superclass::m_DownSampledImage);
  integrator->SetT0(m_T0);
  integrator->SetTN(m_TN);
  integrator->SetJumpTimes(jumpTimes);

  (m_ComputeTrueInverseFlow || m_RegressionFlag) ? integrator->SetComputeTrueInverseFlow()
                                                 : integrator->UnsetComputeTrueInverseFlow();
  m_UseImprovedEuler ? integrator->UseImprovedEuler() : integrator->UseStandardEuler();

  m_UseFastConvolutions ? integrator->SetUseFastConvolutions() : integrator->UnsetUseFastConvolutions();

  integrator->Update();

  m_AdjointPosAt0 = integrator->GetAdjointPosAt(0);
  m_AdjointMomAt0 = integrator->GetAdjointMomAt(0);
  m_AdjointLandmarkPointsAt0 = integrator->GetAdjointLandmarkPointsAt(0);

  delete integrator;
}

template class Diffeos<double,2>;
template class Diffeos<double,3> ;
