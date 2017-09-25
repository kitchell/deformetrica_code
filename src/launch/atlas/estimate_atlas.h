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

/// Core files.
#include "AbstractAtlas.h"
#include "DeterministicAtlas.h"
#include "BayesianAtlas.h"
#include "BayesianAtlasMixture.h"
#include "LdaAtlas.h"

#include "GradientAscent.h"
#include "FastGradientAscent.h"
#include "McmcSaem.h"

#include "SrwMhwgSampler.h"
#include "MalaSampler.h"
#include "AmalaSampler.h"

/// Support files.
#include "ProbabilityDistributions.h"
#include "LinearAlgebra.h"
#include "DeformableObjectReader.h"
#include "SimpleTimer.h"

/// Input-output files.
#include <src/io/XmlDataSet.hpp>
#include "SparseDiffeoParametersXMLFile.h"
#include "DeformableObjectParametersXMLFile.h"

#ifndef DEFORMETRICA_CONFIG
  #include "DeformetricaConfig.h"
#endif

#if ITK_VERSION_MAJOR >= 4
  #include <itkFFTWGlobalConfiguration.h>
#endif

namespace def {
namespace io {
  class XmlModel;
}
};

template<unsigned int Dimension>
void estimateAtlas(std::shared_ptr<const def::io::XmlModel> xml_model);
