/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the MIT License. This file is also distributed     *
*    under the terms of the Inria Non-Commercial License Agreement.                    *
*                                                                                      *
*                                                                                      *
****************************************************************************************/
#ifndef _KernelFactory_h
#define _KernelFactory_h

#include "itkSimpleFastMutexLock.h"
#include "KernelType.h"
#include "ExactKernel.h"

#include <memory>

#include "LinearAlgebra.h"

using namespace def::algebra;

/**
 *	\brief      A kernel factory.
 *
 *	\copyright  Inria and the University of Utah
 *	\version    Deformetrica 2.0
 *
 *	\details    The KernelFactory class enables to instantiate objects whose type is derived from an
 *              abstract type (it is the principle of the factory method pattern). On top of that, this
 *              class implements the singleton pattern, that is to say that you can only instantiate
 *              one factory.
 */

template<class ScalarType, unsigned int PointDim>
class KernelFactory {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Exact kernel type.
  typedef ExactKernel<ScalarType, PointDim> KernelBaseType;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
/*	/// Returns the type of the kernel.
	KernelEnumType GetKernelCode() const { return m_WhichKernel; }
	/// Sets the type of the kernel to ExactKernel.
	void UseExactKernel() { m_Mutex.Lock(); m_WhichKernel = Exact; m_Mutex.Unlock();}
#ifdef USE_CUDA
	/// Sets the type of the kernel to CUDAExactKernel.
	void UseCUDAExactKernel() { m_Mutex.Lock(); m_WhichKernel = CUDAExact; m_Mutex.Unlock();}
#endif
	/// Sets the type of the kernel to P3MKernel.
	void UseP3MKernel() { m_Mutex.Lock(); m_WhichKernel = P3M; m_Mutex.Unlock();}*/
  /// See AbstractDeformations::GetDataDomain() for details.
  static MatrixType GetDataDomain() { return m_DataDomain; }
  /// See AbstractDeformations::SetDataDomain() for details.
  static void SetDataDomain(const MatrixType &DD) {
    m_Mutex.Lock();
    m_DataDomain = DD;
    m_Mutex.Unlock();
  }
  /// See AbstractDeformations::GetWorkingSpacingRatio() for details.
  static ScalarType GetWorkingSpacingRatio() { return m_WorkingSpacingRatio; }
  /// See AbstractDeformations::SetWorkingSpacingRatio() for details.
  static void SetWorkingSpacingRatio(ScalarType d) {
    m_Mutex.Lock();
    m_WorkingSpacingRatio = d;
    m_Mutex.Unlock();
  }
  /// See AbstractDeformations::GetPaddingFactor() for details.
  static ScalarType GetPaddingFactor() { return m_PaddingFactor; }
  /// See AbstractDeformations::SetPaddingFactor() for details.
  static void SetPaddingFactor(ScalarType d) {
    m_Mutex.Lock();
    m_PaddingFactor = d;
    m_Mutex.Unlock();
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Instantiates an object of KernelFactory type with the Singleton strategy.
  static KernelFactory<ScalarType, PointDim> *Instantiate();
  /// Deletes the kernel factory.
  static void Delete();

  /// Returns the instance of the object, NULL in case of error.
  std::shared_ptr<KernelBaseType> CreateKernelObject(KernelEnumType kernelType);

  /// Returns the instance of the object, NULL in case of error.
  std::shared_ptr<KernelBaseType> CreateKernelObject(KernelEnumType kernelType, MatrixType DataDomain);

  /// Returns the instance of the object, NULL in case of error.
  std::shared_ptr<KernelBaseType> CreateKernelObject(
      const MatrixType &X,
      const MatrixType &W,
      ScalarType h, KernelEnumType kernelType);

// NOTE: make sure to assign the target image when using P3M
//   void SetWorkingImage(ImageType* img);
// ImagePointer GetWorkingImage() const { return m_WorkingImage; }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  KernelFactory();
  ~KernelFactory();
// void DetermineP3MGrid(
//   ImagePointType& outOrigin,
//   ImageSpacingType& outSpacing,
//   ImageSizeType& outSize,
//   ImageType* img, ScalarType h);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////
//	/// Type of the kernel.
//	KernelEnumType m_WhichKernel;
  // Information of image domain for heuristically determining parameters
  // ImagePointer m_WorkingImage;

 private:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// See AbstractDeformations::m_DataDomain for details.
  static MatrixType m_DataDomain;
  /// See P3MKernel::m_WorkingSpacingRatio for details.
  static ScalarType m_WorkingSpacingRatio;
  /// See P3MKernel::m_PaddingFactor for details.
  static ScalarType m_PaddingFactor;
  ///	Object used to perform mutex (important for multithreaded programming).
  static itk::SimpleFastMutexLock m_Mutex;
  /// Boolean which enables to instantiate once a kernel factory.
  static bool m_IsInstantiated;
  /// The unique kernel factory.
  static KernelFactory *m_SingletonInstance;
}; /* class KernelFactory */

#endif /* _KernelFactory_h */
