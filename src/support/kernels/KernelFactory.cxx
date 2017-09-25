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

#include "KernelType.h"
#include "KernelFactory.h"
#include "ExactKernel.h"
#include "P3MKernel.h"
#include "Compact.h"

#ifdef USE_CUDA
#include "CUDAExactKernel.h"
#endif

#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization :
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, unsigned int PointDim>
KernelFactory<ScalarType, PointDim> *
KernelFactory<ScalarType, PointDim>::m_SingletonInstance = 0;

template<class ScalarType, unsigned int PointDim>
itk::SimpleFastMutexLock
KernelFactory<ScalarType, PointDim>::m_Mutex;

template<class ScalarType, unsigned int PointDim>
bool
KernelFactory<ScalarType, PointDim>::m_IsInstantiated = false;

template<class ScalarType, unsigned int PointDim>
MatrixType
KernelFactory<ScalarType, PointDim>::m_DataDomain;

template<class ScalarType, unsigned int PointDim>
ScalarType
KernelFactory<ScalarType, PointDim>::m_PaddingFactor = 0;

template<class ScalarType, unsigned int PointDim>
ScalarType
KernelFactory<ScalarType, PointDim>::m_WorkingSpacingRatio = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, unsigned int PointDim>
KernelFactory<ScalarType, PointDim>
::KernelFactory() {}
template<class ScalarType, unsigned int PointDim>
KernelFactory<ScalarType, PointDim>
::~KernelFactory() {}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, unsigned int PointDim>
KernelFactory<ScalarType, PointDim> *
KernelFactory<ScalarType, PointDim>
::Instantiate() {
  if (!m_IsInstantiated) {
    m_Mutex.Lock();
    if (m_SingletonInstance == 0)
      m_SingletonInstance = new KernelFactory<ScalarType, PointDim>();
    m_Mutex.Unlock();
    m_IsInstantiated = true;
  }
  return m_SingletonInstance;
}

template<class ScalarType, unsigned int PointDim>
void
KernelFactory<ScalarType, PointDim>
::Delete() {
  if (m_IsInstantiated) {
    m_Mutex.Lock();
    if (m_SingletonInstance != 0) {
      delete m_SingletonInstance;
      m_SingletonInstance = 0;
    }
    m_Mutex.Unlock();
    m_IsInstantiated = false;
  }
}

template<class ScalarType, unsigned int PointDim>
std::shared_ptr<typename KernelFactory<ScalarType, PointDim>::KernelBaseType>
KernelFactory<ScalarType, PointDim>
::CreateKernelObject(KernelEnumType kernelType) {
  // if (kernelType == null)
  // 	throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");

  switch (kernelType) {
    case Exact: return std::make_shared<ExactKernel<ScalarType, PointDim>>();
#ifdef USE_CUDA
    case CUDAExact:
        return std::make_shared<CUDAExactKernel<ScalarType, PointDim>>();
#endif
    case P3M: {
      typedef P3MKernel<ScalarType, PointDim> P3MKernelType;
      std::shared_ptr<P3MKernelType> obj = std::make_shared<P3MKernelType>();
      if (m_DataDomain.size())
        obj->SetDataDomain(this->GetDataDomain());

      obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
      obj->SetPaddingFactor(this->GetPaddingFactor());
      return obj;
    }
    case COMPACT: return std::make_shared<Compact<ScalarType, PointDim>>();
    default: throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");
  }
}

template<class ScalarType, unsigned int PointDim>
std::shared_ptr<typename KernelFactory<ScalarType, PointDim>::KernelBaseType>
KernelFactory<ScalarType, PointDim>
::CreateKernelObject(KernelEnumType kernelType, MatrixType DataDomain) {
  // if (kernelType == null)
  // 	throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");

  switch (kernelType) {
    case Exact: return std::make_shared<ExactKernel<ScalarType, PointDim>>();
#ifdef USE_CUDA
    case CUDAExact:
        return std::make_shared<CUDAExactKernel<ScalarType, PointDim>>();
#endif
    case P3M: {
      typedef P3MKernel<ScalarType, PointDim> P3MKernelType;
      std::shared_ptr<P3MKernelType> obj = std::make_shared<P3MKernelType>();
      obj->SetDataDomain(DataDomain);
      obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
      obj->SetPaddingFactor(this->GetPaddingFactor());
      return obj;
    }
    case COMPACT: return std::make_shared<Compact<ScalarType, PointDim>>();
    default: throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");
  }
}

template<class ScalarType, unsigned int PointDim>
std::shared_ptr<typename KernelFactory<ScalarType, PointDim>::KernelBaseType>
KernelFactory<ScalarType, PointDim>
::CreateKernelObject(
    const MatrixType &X,
    const MatrixType &W,
    ScalarType h,
    KernelEnumType kernelType) {
  switch (kernelType) {
    case Exact: return std::make_shared<ExactKernel<ScalarType, PointDim>>(X, W, h);
#ifdef USE_CUDA
    case CUDAExact:
        return std::make_shared<CUDAExactKernel<ScalarType, PointDim>>(X, W, h);
#endif
    case P3M: {
      typedef P3MKernel<ScalarType, PointDim> P3MKernelType;
      std::shared_ptr<P3MKernelType> obj = std::make_shared<P3MKernelType>(X, W, h);
      if (m_DataDomain.size())
        obj->SetDataDomain(this->GetDataDomain());
      obj->SetWorkingSpacingRatio(this->GetWorkingSpacingRatio());
      obj->SetPaddingFactor(this->GetPaddingFactor());
      return obj;
    }
    case COMPACT: return std::make_shared<Compact<ScalarType, PointDim>>();
    default: throw std::runtime_error("In KernelFactory::CreateKernelObject() - The type of the kernel is unknown");
  }
}

template class KernelFactory<double, 2>;
template class KernelFactory<double, 3>;
