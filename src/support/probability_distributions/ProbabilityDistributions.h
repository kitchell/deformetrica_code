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

#include "NormalDistribution.h"
#include "AbstractNormalDistribution.h"
#include "NormalDistribution.h"
#include "MultiScalarNormalDistribution.h"
#include "DisplacementFieldNormalDistribution.h"
#include "UniformDistribution.h"
#include "NormalDistribution.h"
#include "InverseWishartDistribution.h"
#include "MultiScalarInverseWishartDistribution.h"
#include "DirichletDistribution.h"
#include "AutomaticRelevanceDeterminationDistribution.h"

namespace def {
namespace proba {

#ifdef USE_DOUBLE_PRECISION
typedef double ScalarType;
#else
typedef float ScalarType;
#endif

typedef AbstractProbabilityDistribution<ScalarType> ProbabilityDistributionType;
typedef std::map<std::string, std::shared_ptr<ProbabilityDistributionType>> ProbabilityDistributionMapType;
typedef AbstractNormalDistribution<ScalarType> AbstractNormalDistributionType;
typedef std::map<std::string, std::shared_ptr<AbstractNormalDistributionType>> AbstractNormalDistributionMapType;
typedef NormalDistribution<ScalarType> NormalDistributionType;
typedef std::map<std::string, std::shared_ptr<NormalDistributionType>> NormalDistributionMapType;
typedef MultiScalarNormalDistribution<ScalarType> MultiScalarNormalDistributionType;
typedef std::map<std::string, std::shared_ptr<MultiScalarNormalDistributionType>> MultiScalarNormalDistributionMapType;
typedef DisplacementFieldNormalDistribution<ScalarType> DisplacementFieldNormalDistributionType;
typedef std::map<std::string, std::shared_ptr<DisplacementFieldNormalDistributionType>> DisplacementFieldNormalDistributionMapType;
typedef UniformDistribution<ScalarType> UniformDistributionType;
typedef InverseWishartDistribution<ScalarType> InverseWishartDistributionType;
typedef MultiScalarInverseWishartDistribution<ScalarType> MultiScalarInverseWishartDistributionType;
typedef DirichletDistribution<ScalarType> DirichletDistributionType;
typedef AutomaticRelevanceDeterminationDistribution<ScalarType> AutomaticRelevanceDeterminationDistributionType;
}
}
