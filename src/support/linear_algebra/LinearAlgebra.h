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

#include "ArmadilloVectorWrapper.h"
#include "ArmadilloMatrixWrapper.h"
#include "MatrixListWrapper.h"
#include "LinearVariableWrapper.h"
#include "LinearVariablesWrapper.h"
#include "LinearVariableMapWrapper.h"
#include "LinearVariablesMapWrapper.h"

namespace def {
namespace algebra {

#ifdef USE_DOUBLE_PRECISION
typedef double ScalarType;
#else
typedef float ScalarType;
#endif

typedef ArmadilloVectorWrapper<ScalarType> VectorType;
typedef ArmadilloMatrixWrapper<ScalarType> MatrixType;

typedef MatrixListWrapper<ScalarType> MatrixListType;

typedef LinearVariableWrapper<ScalarType> LinearVariableType;
typedef LinearVariablesWrapper<ScalarType> LinearVariablesType;

typedef LinearVariableMapWrapper<ScalarType> LinearVariableMapType;
typedef LinearVariablesMapWrapper<ScalarType> LinearVariablesMapType;

}
}
