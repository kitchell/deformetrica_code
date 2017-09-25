#pragma once

namespace def {
namespace io {
class XmlModel;
}
};

#include "SparseDiffeoParameters.h"
#include "DeformableObjectParameters.h"

#include "LinearAlgebra.h"
using namespace def::algebra;

template <class ScalarType, unsigned int Dimension>
void deform(SparseDiffeoParameters* paramDiffeos,
            bool useInverseFlow,
            const char* CP_fn,
            MatrixType MOM0_i,
            int numObjects,
            std::vector<DeformableObjectParameters::Pointer> paramObjectsList,
            const std::vector<std::string> objectfnList);

template <unsigned int Dimension>
void deformation(std::shared_ptr<const def::io::XmlModel> xml_model);