/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/
#include "itkTranslationTransform.h"


#include "MutualInformationImage.h"
#include <src/core/observations/deformable_objects/geometries/AbstractGeometry.h>

template<class ScalarType, unsigned int Dimension>
MutualInformationImage<ScalarType, Dimension>
::MutualInformationImage() : Superclass() {
  this->SetMutualInformationImageType();
}

template<class ScalarType, unsigned int Dimension>
MutualInformationImage<ScalarType, Dimension>
::MutualInformationImage(const MutualInformationImage<ScalarType, Dimension> &o) : Superclass(o) {}

template<class ScalarType, unsigned int Dimension>
MutualInformationImage<ScalarType, Dimension>
::MutualInformationImage(const MutualInformationImage<ScalarType, Dimension> &ex, const MatrixType &IP) : Superclass(ex, IP) {
  this->SetMutualInformationImageType();
  this->Update();
}

template<class ScalarType, unsigned int Dimension>
MutualInformationImage<ScalarType, Dimension>
::~MutualInformationImage() {};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, unsigned int Dimension>
void
MutualInformationImage<ScalarType, Dimension>
::Update() {
  if (this->IsModified()) {
    Superclass::Update();
  }
  this->UnSetModified();
}

template<class ScalarType, unsigned int Dimension>
ScalarType MutualInformationImage<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<AbstractGeometryType> target)
{
  if (this->GetType() != target->GetType())
    throw std::runtime_error("Abstract Geometries types mismatched");

  const std::shared_ptr<const MutualInformationImage> targ = std::static_pointer_cast<const MutualInformationImage>(target);

  typedef itk::TranslationTransform<ScalarType, Dimension> TranslationTransformType; // This cannot be float for some reason?
  typename TranslationTransformType::Pointer transform = TranslationTransformType::New();

  typedef itk::MutualInformationImageToImageMetric<ImageType,ImageType> MetricType;

  typename MetricType::Pointer metric = MetricType::New();

  metric->SetTransform(transform);
  metric->SetFixedImageStandardDeviation(0.3);
  metric->SetMovingImageStandardDeviation(0.3);
  metric->SetFixedImage(this->GetImage());
  metric->SetMovingImage(targ->GetImage());
  metric->SetFixedImageRegion(this->GetImage()->GetLargestPossibleRegion());

  typename itk::LinearInterpolateImageFunction<ImageType, ScalarType>::Pointer interpolator = itk::LinearInterpolateImageFunction<ImageType, double>::New();
  interpolator->SetInputImage(this->GetImage());
  metric->SetInterpolator(interpolator);

  typename TranslationTransformType::ParametersType parameters;
  parameters.SetSize(2);
  parameters.Fill(0);

  ScalarType match = metric->GetValue(parameters);

  return match;

}

// compute 2 (I_0(y(0)) - I_1) nabla_{y(0)}I_0
template<class ScalarType, unsigned int Dimension>
MatrixType
MutualInformationImage<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<AbstractGeometryType> target) {
  throw std::runtime_error("ComputeMatchGradient not implemented for MutualInformationImage");
}

template class MutualInformationImage<double,2>;
template class MutualInformationImage<double,3>;