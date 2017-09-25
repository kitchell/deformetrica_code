/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "SparseDiffeoParameters.h"
#include "itksys/SystemTools.hxx"


////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

SparseDiffeoParameters
::SparseDiffeoParameters() {
  m_KernelWidth = 0.0;

  m_KernelType = "cudaexact";
  // 1/5 gives a relative approximation error of about 5%, use 0.3 for increased speed
  m_P3MWorkingSpacingRatio = 0.2;
  // enlarge grids by 3 x kernelwidth to avoid side effects (FFTs have circular boundary conditions). It is also used to define a bounding box
  m_P3MPaddingFactor = 3.0f;

  m_T0 = 0.0;
  m_TN = 1.0;
  m_NumberOfTimePoints = 10;
  m_InitialCPSpacing = 0.0;

  ///For abc_sampling:
  m_NumberOfCps = 1;
  m_NumberOfSamples = 10;
  m_MomentaPropositionVariance = 1.;
  m_TrainingSetSize = 10;

  ///For parallel transport
  m_MatchingNumberOfTimePoints = 10;
  m_UseExpParallelization = true;
  m_MatchingT0 = 0.;
  m_MatchingTN = 1.;
  m_InitialCPPositionForTransport_fn = "";
  m_InitialMomentaForTransport_fn = "";

  m_OptimizationMethodType = "FastGradientAscent";

  m_CovarianceMomenta_Normalized_Hyperparameter = 0.2;
  m_CovarianceMomenta_Prior_Inverse_fn = "";
  m_CovarianceMomentaInverse_fn = "";
  m_ModelType = "Undefined";
  m_MaximumNumberOfClasses = 1;

  m_ModelName = "";

  m_UseTempering = false;
  m_InitialTemperature = 100.0;
  m_TemperingDurationRatio = 0.5;

  m_FreezeCP = false;
  m_FreezeTemplate = false;

  m_InitialCPPosition_fn = "";
  m_InitialMomenta_fn = "";
  m_ModulationMatrix_fn = "";

  m_ReferenceTime = sqrt(- 1.0);
  m_TimeShiftVariance = - 1.0;
  m_LogAccelerationVariance = - 1.0;

  m_SmoothingKernelWidthRatio = 1;

  m_MaxIterations = 100;
  m_MaxLineSearchIterations = 10;

  m_StepExpand = 1.2;
  m_StepShrink = 0.5;

  m_AdaptiveTolerance = 1e-4;
  m_InitialStepSize = 0.001;

  m_ComputeTrueInverseFlow = Default;
  m_UseImplicitEuler = false;
  m_UseImplicitEuler = true;
  m_OptimizeInitialControlPoints = false;
  m_UseFastConvolutions = false;
  m_MultivariateLineSearch = true;
  m_WriteFullTrajectories = false;

  m_NumberOfThreads = 1;

  m_ConcentrationOfTimePointsForReferenceGeodesic = 10;
  m_NumberOfTimePointsForExponentiation = 10;
  m_MarginOnGeodesicLength = 1.25;

  m_NumberOfSources = 4;

  m_LogAccelerationRandomEffectStd = -1;
  m_ReferenceTimePriorMean = sqrt(-1);
  m_ReferenceTimePriorStd = -1;

  m_PrintEveryNIters = 1;
  m_SaveEveryNIters = 1;

  m_PrintAcceptanceRatesWindow = 100;
  m_AdaptiveAcceptanceRatesWindow = 10;
  m_AdaptiveAcceptanceRatesTarget = 30;

  m_SrwProposalBlocksize = 0;

  m_SrwTemplateDataProposalStd = 0.01;
  m_SrwTemplateDataProposalKernelWidth = -1;
  m_SrwTemplateDataProposalNumberOfControlPoints = 10;

  m_SrwControlPointsProposalStd = 0.01;
  m_SrwMomentaProposalStd = 0.01;
  m_SrwModulationMatrixProposalStd = 0.01;
  m_SrwSourcesProposalStd = 0.01;
  m_SrwLogAccelerationProposalStd = 0.01;
  m_SrwTimeShiftProposalStd = 0.01;
  m_SrwMixtureCoefficientsProposalStd = 0.01;

  m_MalaControlPointsScale = 0.001;
  m_MalaMomentaScale = 0.001;

  m_AmalaControlPointsMeanScale = 0.001;
  m_AmalaControlPointsCovarianceScale = 0.002;
  m_AmalaMomentaMeanScale = 0.001;
  m_AmalaMomentaCovarianceScale = 0.002;
  m_AmalaControlPointsRegularization = 0.0001;
  m_AmalaMomentaRegularization = 0.0001;

  ///For LDA Model:
  m_IntraClassPCADimension = 5;
  m_InitialG_fn = "";
  m_InitialBeta_fn = "";
  m_InitialF_fn = "";
  m_InitialAlpha_fn = "";
  m_OnlyWrite = false;

}

SparseDiffeoParameters
::~SparseDiffeoParameters() {

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

void
SparseDiffeoParameters
::Update() {
  if (m_KernelWidth < 1e-20) {
    std::cout << "Warning : no value set for Kernel Width. Defaulted to 1." << std::endl;
    m_KernelWidth = 1;
  }

  if (!(m_InitialCPPosition_fn.size()) && m_InitialCPSpacing < 1e-20) {
    std::cout << "No initial CP spacing given: using diffeo kernel width of " << m_KernelWidth
              << "." << std::endl;
    m_InitialCPSpacing = m_KernelWidth;
  }
}

bool
SparseDiffeoParameters
::CheckValues() const {
  if (m_KernelWidth < 1e-20)
    return false;

  if (m_NumberOfTimePoints < 2)
    return false;

  if (m_MaxIterations < 1)
    return false;
  if (m_MaxLineSearchIterations < 1)
    return false;

  if (m_AdaptiveTolerance <= 0.0)
    return false;
  if (m_InitialStepSize < 1e-20)
    return false;

  if (m_InitialCPSpacing < 1e-20)
    return false;

  //  ********************** CHANGE *******************
  if (m_CovarianceMomenta_Normalized_Hyperparameter < 0.0)
    return false;
  //  ********************** CHANGE *******************

  return true;
}

void
SparseDiffeoParameters
::PrintSelf(std::ostream &os) {
  os << "Kernel width = " << m_KernelWidth << std::endl;
  os << std::endl;
  os << "Kernel type = " << m_KernelType << std::endl;
  os << std::endl;
  os << "T0 = " << m_T0 << "   Tn = " << m_TN << std::endl;
  os << "Number of time points = " << m_NumberOfTimePoints << std::endl;
  os << "Initial CP spacing = " << m_InitialCPSpacing << std::endl;
  os << "Freeze CP = " << (m_FreezeCP ? "On" : "Off") << std::endl;
  os << "Initial set of CP loaded from " << m_InitialCPPosition_fn << std::endl;
  os << "Initial momenta loaded from " << m_InitialMomenta_fn << std::endl;
  os << std::endl;
  os << "Freeze Template = " << (m_FreezeTemplate ? "On" : "Off") << std::endl;
  os << std::endl;
  os << std::endl;
  os << "SmoothingKernelWidthRatio = " << m_SmoothingKernelWidthRatio << std::endl;
  os << "Optimization method: " << m_OptimizationMethodType << std::endl;
  os << "Covariance Momenta inverse loaded from " << m_CovarianceMomentaInverse_fn << std::endl;
  os << "Using Fast convolutions :" << m_UseFastConvolutions << std::endl;
  os << "Using improved euler " << m_UseImprovedEuler << std::endl;
  os << "Covariance Momenta Normalized Hyperparameter: " << m_CovarianceMomenta_Normalized_Hyperparameter << std::endl;
//	os << "Bayesian Framework: " << (m_BayesianFramework?"On":"Off") << std::endl;
  os << "Model type = " << m_ModelType << std::endl;
  os << "Max number of class = " << m_MaximumNumberOfClasses << std::endl;
  os << "Inverse of the Prior for Covariance Momenta loaded from " << m_CovarianceMomenta_Prior_Inverse_fn << std::endl;
  os << "Intra class PCA Dimension " << m_IntraClassPCADimension << std::endl;
  os << "Max descent iterations = " << m_MaxIterations << std::endl;
  os << "Max line search iterations = " << m_MaxLineSearchIterations << std::endl;
  os << "Step expand = " << m_StepExpand << std::endl;
  os << "Step shrink = " << m_StepShrink << std::endl;
  os << "Adaptive tolerance = " << m_AdaptiveTolerance << std::endl;
  os << "Initial step size = " << m_InitialStepSize << std::endl;
  os << "Number of threads = " << m_NumberOfThreads << std::endl;
  os << "Compute True Inverse Flow (for images) = " << (m_ComputeTrueInverseFlow ? "On" : "Off") << std::endl;
  os << "P3M working spacing ratio (for kernels of P3M type) = " << m_P3MWorkingSpacingRatio << std::endl;
  os << "P3M padding factor (for kernels of P3M type) = " << m_P3MPaddingFactor << std::endl;
  os << std::endl;
}

