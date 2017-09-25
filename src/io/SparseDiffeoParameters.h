/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#pragma once

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>

#include <limits>
#include <cstddef>

class SparseDiffeoParameters : public itk::Object {

 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Typedefs :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef SparseDiffeoParameters Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef enum {
    Default = 0,
    On = 1,
    Off = 2
  } BooleanOptionType;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  SparseDiffeoParameters();
  ~SparseDiffeoParameters();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  // Compulsory parameters
  itkGetMacro(KernelWidth, double);
  itkSetMacro(KernelWidth, double);

  // Optional parameters...
  // ... for the diffeos
  itkGetMacro(KernelType, std::string);
  itkSetMacro(KernelType, std::string);

  itkGetMacro(NumberOfTimePoints, unsigned
      int);
  itkSetMacro(NumberOfTimePoints, unsigned
      int);

  itkGetMacro(NumberOfCps, unsigned int);
  itkSetMacro(NumberOfCps, unsigned int);

  itkGetMacro(NumberOfSamples, unsigned int);
  itkSetMacro(NumberOfSamples, unsigned int);


  itkGetMacro(MomentaPropositionVariance, double);
  itkSetMacro(MomentaPropositionVariance, double);

  itkGetMacro(TrainingSetSize, unsigned int);
  itkSetMacro(TrainingSetSize, unsigned int);

  itkGetMacro(T0, double);
  itkSetMacro(T0, double);
  itkGetMacro(TN, double);
  itkSetMacro(TN, double);

  itkGetMacro(InitialCPSpacing, double);
  itkSetMacro(InitialCPSpacing, double);

  itkGetMacro(P3MWorkingSpacingRatio, double);
  itkSetMacro(P3MWorkingSpacingRatio, double);

  itkGetMacro(P3MPaddingFactor, double);
  itkSetMacro(P3MPaddingFactor, double);


//	// ... for the gradient ascent
//	void SetUseFISTA() { m_UseFISTA = true; }
//	void UnsetUseFISTA() { m_UseFISTA = false; }
//	bool UseFISTA() { return m_UseFISTA; }

  // ... for the optimization method
  itkGetMacro(OptimizationMethodType, std::string);
  itkSetMacro(OptimizationMethodType, std::string);

//	void SetAdaptiveDataSigma() { m_AdaptiveDataSigma = true; }
//	void UnsetAdaptiveDataSigma() { m_AdaptiveDataSigma = false; }
//	bool AdaptiveDataSigma() { return m_AdaptiveDataSigma; }

  //  ********************** CHANGE *******************
  void SetCovarianceMomenta_Normalized_Hyperparameter(double d) {
    m_CovarianceMomenta_Normalized_Hyperparameter = d;
  }
  double GetCovarianceMomenta_Normalized_Hyperparameter() { return m_CovarianceMomenta_Normalized_Hyperparameter; }

//	void SetBayesianFramework() { m_BayesianFramework = true; }
//	void UnsetBayesianFramework() { m_BayesianFramework = false; }
//	bool BayesianFramework() { return m_BayesianFramework; }
  itkGetMacro(ModelType, std::string);
  itkSetMacro(ModelType, std::string);

  itkGetMacro(MaximumNumberOfClasses, int);
  itkSetMacro(MaximumNumberOfClasses, int);

  itkGetMacro(CovarianceMomentaInverse_fn, std::string);
  itkSetMacro(CovarianceMomentaInverse_fn, std::string);

  itkGetMacro(CovarianceMomenta_Prior_Inverse_fn, std::string);
  itkSetMacro(CovarianceMomenta_Prior_Inverse_fn, std::string);

  itkGetMacro(ModelName, std::string);
  itkSetMacro(ModelName, std::string);

  itkGetMacro(InitialCPPosition_fn, std::string);
  itkSetMacro(InitialCPPosition_fn, std::string);

  itkGetMacro(InitialMomenta_fn, std::string);
  itkSetMacro(InitialMomenta_fn, std::string);

  itkGetMacro(ModulationMatrix_fn, std::string);
  itkSetMacro(ModulationMatrix_fn, std::string);

  itkGetMacro(ReferenceTime, double);
  itkSetMacro(ReferenceTime, double);

  itkGetMacro(TimeShiftVariance, double);
  itkSetMacro(TimeShiftVariance, double);

  itkGetMacro(LogAccelerationVariance, double);
  itkSetMacro(LogAccelerationVariance, double);

  itkGetMacro(InitialCPPositionForTransport_fn, std::string);
  itkSetMacro(InitialCPPositionForTransport_fn, std::string);

  itkGetMacro(InitialMomentaForTransport_fn, std::string);
  itkSetMacro(InitialMomentaForTransport_fn, std::string);

  itkGetMacro(MatchingNumberOfTimePoints, unsigned int);
  itkSetMacro(MatchingNumberOfTimePoints, unsigned int);

  itkGetMacro(MatchingT0, double);
  itkSetMacro(MatchingT0, double);

  itkGetMacro(MatchingTN, double);
  itkSetMacro(MatchingTN, double);

  void SetUseExpParallelization() {m_UseExpParallelization = true;}
  void UnsetUseExpParallelization() {m_UseExpParallelization = false;}
  bool UseExpParallelization() {return m_UseExpParallelization;}

  void SetUseTempering() {m_UseTempering = true;}
  void UnsetUseTempering() {m_UseTempering = false;}
  bool UseTempering() {return m_UseTempering;}

  itkGetMacro(InitialTemperature, double);
  itkSetMacro(InitialTemperature, double);

  itkGetMacro(TemperingDurationRatio, double);
  itkSetMacro(TemperingDurationRatio, double);

  void SetFreezeCP() {m_FreezeCP = true;}
  void UnsetFreezeCP() {m_FreezeCP = false;}
  bool FreezeCP() {return m_FreezeCP;}

  void SetFreezeTemplate() { m_FreezeTemplate = true; }
  void UnsetFreezeTemplate() { m_FreezeTemplate = false; }
  bool FreezeTemplate() { return m_FreezeTemplate; }

  itkGetMacro(SmoothingKernelWidthRatio, double);
  itkSetMacro(SmoothingKernelWidthRatio, double);

  itkGetMacro(MaxIterations, unsigned
      int);
  itkSetMacro(MaxIterations, unsigned
      int);

  itkGetMacro(MaxLineSearchIterations, unsigned
      int);
  itkSetMacro(MaxLineSearchIterations, unsigned
      int);

  itkGetMacro(StepExpand, double);
  itkSetMacro(StepExpand, double);

  itkGetMacro(StepShrink, double);
  itkSetMacro(StepShrink, double);

  itkGetMacro(AdaptiveTolerance, double);
  itkSetMacro(AdaptiveTolerance, double);

  itkGetMacro(InitialStepSize, double);
  itkSetMacro(InitialStepSize, double);

  itkGetMacro(NumberOfThreads, unsigned
      int);
  itkSetMacro(NumberOfThreads, unsigned
      int);

  itkGetMacro(ConcentrationOfTimePointsForReferenceGeodesic, unsigned
      int);
  itkSetMacro(ConcentrationOfTimePointsForReferenceGeodesic, unsigned
      int);
  itkGetMacro(NumberOfTimePointsForExponentiation, unsigned
      int);
  itkSetMacro(NumberOfTimePointsForExponentiation, unsigned
      int);
  itkGetMacro(MarginOnGeodesicLength, double);
  itkSetMacro(MarginOnGeodesicLength, double);

  itkGetMacro(NumberOfSources, unsigned int);
  itkSetMacro(NumberOfSources, unsigned int);

  itkGetMacro(LogAccelerationRandomEffectStd, double);
  itkSetMacro(LogAccelerationRandomEffectStd, double);

  itkGetMacro(ReferenceTimePriorMean, double);
  itkSetMacro(ReferenceTimePriorMean, double);
  itkGetMacro(ReferenceTimePriorStd, double);
  itkSetMacro(ReferenceTimePriorStd, double);

  void SetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = On; }
  void UnsetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = Off; }
  BooleanOptionType ComputeTrueInverseFlow() { return m_ComputeTrueInverseFlow; }

  void SetUseImplicitEuler() { m_UseImplicitEuler = true; }
  void UnsetUseImplicitEuler() { m_UseImplicitEuler = false; }
  bool UseImplicitEuler() { return m_UseImplicitEuler; }

  void SetUseImprovedEuler() { m_UseImprovedEuler = true; }
  void UnsetUseImprovedEuler() { m_UseImprovedEuler = false; }
  bool UseImprovedEuler() { return m_UseImprovedEuler; }

  void SetOptimizeInitialControlPoints() { m_OptimizeInitialControlPoints = true; }
  void UnsetOptimizeInitialControlPoints() { m_OptimizeInitialControlPoints = false; }
  bool OptimizeInitialControlPoints() { return m_OptimizeInitialControlPoints; }

  void SetUseFastConvolutions() { m_UseFastConvolutions = true; }
  void UnsetUseFastConvolutions() { m_UseFastConvolutions = false; }
  bool UseFastConvolutions() { return m_UseFastConvolutions; }

  void SetMultivariateLineSearch() { m_MultivariateLineSearch = true; }
  void UnsetMultivariateLineSearch() { m_MultivariateLineSearch = false; }
  bool MultivariateLineSearch() { return m_MultivariateLineSearch; }

  void SetWriteFullTrajectories() { m_WriteFullTrajectories = true; }
  void UnsetWriteFullTrajectories() { m_WriteFullTrajectories = false; }
  bool WriteFullTrajectories() { return m_WriteFullTrajectories; }

  void SetUseForwardShooting() { m_UseForwardShooting = true; }
  void UnsetUseForwardShooting() { m_UseForwardShooting = false; }
  bool UseForwardShooting() { return m_UseForwardShooting; }

  itkGetMacro(PrintEveryNIters, unsigned
      int);
  itkSetMacro(PrintEveryNIters, unsigned
      int);

  itkGetMacro(SaveEveryNIters, unsigned
      int);
  itkSetMacro(SaveEveryNIters, unsigned
      int);

  itkGetMacro(PrintAcceptanceRatesWindow, unsigned
      int);
  itkSetMacro(PrintAcceptanceRatesWindow, unsigned
      int);

  itkGetMacro(AdaptiveAcceptanceRatesWindow, unsigned
      int);
  itkSetMacro(AdaptiveAcceptanceRatesWindow, unsigned
      int);

  itkGetMacro(AdaptiveAcceptanceRatesTarget, double);
  itkSetMacro(AdaptiveAcceptanceRatesTarget, double);

  itkGetMacro(SrwProposalBlocksize, unsigned
      int);
  itkSetMacro(SrwProposalBlocksize, unsigned
      int);

  itkGetMacro(SrwTemplateDataProposalStd, double);
  itkSetMacro(SrwTemplateDataProposalStd, double);

  itkGetMacro(SrwTemplateDataProposalKernelWidth, double);
  itkSetMacro(SrwTemplateDataProposalKernelWidth, double);

  itkGetMacro(SrwTemplateDataProposalNumberOfControlPoints, unsigned int);
  itkSetMacro(SrwTemplateDataProposalNumberOfControlPoints, unsigned int);

  itkGetMacro(SrwControlPointsProposalStd, double);
  itkSetMacro(SrwControlPointsProposalStd, double);

  itkGetMacro(SrwMomentaProposalStd, double);
  itkSetMacro(SrwMomentaProposalStd, double);

  itkGetMacro(SrwModulationMatrixProposalStd, double);
  itkSetMacro(SrwModulationMatrixProposalStd, double);

  itkGetMacro(SrwSourcesProposalStd, double);
  itkSetMacro(SrwSourcesProposalStd, double);

  itkGetMacro(SrwLogAccelerationProposalStd, double);
  itkSetMacro(SrwLogAccelerationProposalStd, double);

  itkGetMacro(SrwTimeShiftProposalStd, double);
  itkSetMacro(SrwTimeShiftProposalStd, double);

  itkGetMacro(SrwMixtureCoefficientsProposalStd, double);
  itkSetMacro(SrwMixtureCoefficientsProposalStd, double);

  itkGetMacro(MalaControlPointsScale, double);
  itkSetMacro(MalaControlPointsScale, double);

  itkGetMacro(MalaMomentaScale, double);
  itkSetMacro(MalaMomentaScale, double);

  itkGetMacro(AmalaControlPointsMeanScale, double);
  itkSetMacro(AmalaControlPointsMeanScale, double);

  itkGetMacro(AmalaControlPointsCovarianceScale, double);
  itkSetMacro(AmalaControlPointsCovarianceScale, double);

  itkGetMacro(AmalaMomentaMeanScale, double);
  itkSetMacro(AmalaMomentaMeanScale, double);

  itkGetMacro(AmalaMomentaCovarianceScale, double);
  itkSetMacro(AmalaMomentaCovarianceScale, double);

  itkGetMacro(AmalaControlPointsRegularization, double);
  itkSetMacro(AmalaControlPointsRegularization, double);

  itkGetMacro(AmalaMomentaRegularization, double);
  itkSetMacro(AmalaMomentaRegularization, double);

  ///For LDA model
  itkGetMacro(IntraClassPCADimension, unsigned int);
  itkSetMacro(IntraClassPCADimension, unsigned int);

  itkGetMacro(InitialG_fn, std::string);
  itkSetMacro(InitialG_fn, std::string);

  itkGetMacro(InitialBeta_fn, std::string);
  itkSetMacro(InitialBeta_fn, std::string);

  itkGetMacro(InitialAlpha_fn, std::string);
  itkSetMacro(InitialAlpha_fn, std::string);

  itkGetMacro(InitialF_fn, std::string);
  itkSetMacro(InitialF_fn, std::string);

  void SetOnlyWrite() { m_OnlyWrite = true; }
  void UnsetOnlyWrite() { m_OnlyWrite = false; }
  bool OnlyWrite() { return m_OnlyWrite; }

  void SetFreezeF() { m_FreezeF = true; }
  void UnsetFreezeF() { m_FreezeF = false; }
  bool FreezeF() { return m_FreezeF; }

  void SetFreezeG() { m_FreezeG = true; }
  void UnsetFreezeG() { m_FreezeG = false; }
  bool FreezeG() { return m_FreezeG; }


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other public method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  void Update();

  itkNewMacro(Self);

  // Make sure all values are OK
  virtual bool CheckValues() const;

  void PrintSelf(std::ostream &os);

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  double m_KernelWidth;

  std::string m_KernelType;
  double m_T0;
  double m_TN;

  unsigned int m_NumberOfCps;
  unsigned int m_NumberOfSamples;
  double m_MomentaPropositionVariance;
  unsigned int m_TrainingSetSize;

  unsigned int m_NumberOfTimePoints;
  double m_InitialCPSpacing;
  double m_P3MWorkingSpacingRatio;
  double m_P3MPaddingFactor;

  BooleanOptionType m_ComputeTrueInverseFlow;
  bool m_UseImplicitEuler;
  bool m_UseImprovedEuler;
  bool m_OptimizeInitialControlPoints;
  bool m_UseFastConvolutions;
  bool m_MultivariateLineSearch;
  bool m_WriteFullTrajectories;

  //the old +1 or -1 command line input for shootAndFlow
  bool m_UseForwardShooting;

  std::string m_OptimizationMethodType;
//	bool m_AdaptiveDataSigma;

  double m_CovarianceMomenta_Normalized_Hyperparameter;
//	bool m_BayesianFramework;
  std::string m_ModelType;
  int m_MaximumNumberOfClasses;
  std::string m_CovarianceMomenta_Prior_Inverse_fn; // Inverse of the Prior of the Cov Momenta in a Bayesian Framework
  std::string m_CovarianceMomentaInverse_fn; // Matrix to use in place of the kernel in a Deterministic Framework

  std::string m_ModelName;

  bool m_UseTempering;
  double m_InitialTemperature;
  double m_TemperingDurationRatio;

  bool m_FreezeCP;
  bool m_FreezeTemplate;
  std::string m_InitialCPPosition_fn;
  std::string m_InitialMomenta_fn;
  std::string m_ModulationMatrix_fn;

  double m_ReferenceTime;
  double m_TimeShiftVariance;
  double m_LogAccelerationVariance;

  ///For parallel transport
  std::string m_InitialCPPositionForTransport_fn;
  std::string m_InitialMomentaForTransport_fn;
  unsigned int m_MatchingNumberOfTimePoints;
  bool m_UseExpParallelization;
  double m_MatchingT0;
  double m_MatchingTN;


  double m_SmoothingKernelWidthRatio;

  unsigned int m_MaxIterations;
  unsigned int m_MaxLineSearchIterations;

  double m_StepExpand;
  double m_StepShrink;
  double m_AdaptiveTolerance;
  double m_InitialStepSize;

  unsigned int m_NumberOfThreads;

  unsigned int m_ConcentrationOfTimePointsForReferenceGeodesic;
  unsigned int m_NumberOfTimePointsForExponentiation;
  double m_MarginOnGeodesicLength;

  unsigned int m_NumberOfSources;

  double m_LogAccelerationRandomEffectStd;
  double m_ReferenceTimePriorMean;
  double m_ReferenceTimePriorStd;

  unsigned int m_PrintEveryNIters;
  unsigned int m_SaveEveryNIters;

  unsigned int m_PrintAcceptanceRatesWindow;
  unsigned int m_AdaptiveAcceptanceRatesWindow;
  double m_AdaptiveAcceptanceRatesTarget;

  unsigned int m_SrwProposalBlocksize;

  double m_SrwTemplateDataProposalStd;
  double m_SrwTemplateDataProposalKernelWidth;
  unsigned int m_SrwTemplateDataProposalNumberOfControlPoints;

  double m_SrwControlPointsProposalStd;
  double m_SrwMomentaProposalStd;
  double m_SrwModulationMatrixProposalStd;
  double m_SrwSourcesProposalStd;
  double m_SrwTimeShiftProposalStd;
  double m_SrwLogAccelerationProposalStd;
  double m_SrwMixtureCoefficientsProposalStd;

  double m_MalaControlPointsScale;
  double m_MalaMomentaScale;

  double m_AmalaControlPointsMeanScale;
  double m_AmalaControlPointsCovarianceScale;
  double m_AmalaMomentaMeanScale;
  double m_AmalaMomentaCovarianceScale;
  double m_AmalaControlPointsRegularization;
  double m_AmalaMomentaRegularization;

  ///For LDA model
  unsigned int m_IntraClassPCADimension;
  std::string m_InitialG_fn;
  std::string m_InitialBeta_fn;
  std::string m_InitialF_fn;
  std::string m_InitialAlpha_fn;
  bool m_OnlyWrite;
  bool m_FreezeF;
  bool m_FreezeG;
}; /* class SparseDiffeoParameters */

