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

#include <fstream>
#include <sstream>
#include <string>

#include "LinearAlgebra.h"

using namespace def::algebra;

template <class ScalarType>
MatrixType readMatrixDLM(const char* fn);

template <class ScalarType>
std::vector<MatrixType> readMultipleMatrixDLM(const char* fn);

template <class ScalarType>
void writeMatrixDLM(std::string fn, const MatrixType & M);

template <class ScalarType>
void writeMultipleMatrixDLM(std::string fn, MatrixListType const& M);

template<class ScalarType>
void printMatrix(std::string const name, MatrixType const& M);

