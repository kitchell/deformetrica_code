/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#include "MatrixDLM.h"
#include "ParametricImage.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
ParametricImage<ScalarType, Dimension>
::ParametricImage() : Superclass()
{
  this->SetParametricImageType();
}



template <class ScalarType, unsigned int Dimension>
ParametricImage<ScalarType, Dimension>
::ParametricImage(const ParametricImage& other) : Superclass(other)
{
  m_PhotometricControlPoints = other.m_PhotometricControlPoints;
  m_PhotometricWeights = other.m_PhotometricWeights;

  m_PhotometricKernelType = other.m_PhotometricKernelType;
  m_PhotometricKernelWidth = other.m_PhotometricKernelWidth;
  m_PhotometricKernelField = other.m_PhotometricKernelField;

  m_PhotometricCPSpacing = other.m_PhotometricCPSpacing;
  m_NbPhotoCPs = other.m_NbPhotoCPs;
  m_NbVoxels = other.m_NbVoxels;

  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

  if (!other.m_Image.IsNull()) {
    duplicator->SetInputImage(other.m_Image);
    duplicator->Update();
    m_Image = duplicator->GetOutput(); }

  if (!other.m_DownSampledImage.IsNull()) {
    duplicator->SetInputImage(other.m_DownSampledImage);
    duplicator->Update();
    m_DownSampledImage = duplicator->GetOutput(); }

  m_ImageSpatialGradient.resize(Dimension);
  for (unsigned int dim = 0; dim < Dimension; dim++)
  {
    if (!other.m_ImageSpatialGradient[dim].IsNull()) {
      duplicator->SetInputImage(other.m_ImageSpatialGradient[dim]);
      duplicator->Update();
      m_ImageSpatialGradient[dim] = duplicator->GetOutput(); }
  }

  m_LIImage = other.m_LIImage;

  m_PermutationAxes = other.m_PermutationAxes;
  m_FlipAxes = other.m_FlipAxes;

  this->SetParametricImageType();
  Update();
}


template <class ScalarType, unsigned int Dimension>
ParametricImage<ScalarType, Dimension>
::ParametricImage(const ParametricImage& example, const MatrixType& downSampledImageMap) : Superclass(example, downSampledImageMap)
{
  if (example.m_Image.IsNull() || example.m_DownSampledImage.IsNull())
    throw std::runtime_error("ParametricImage : Cannot re-sample image if no image has been set.");

  /// Copy example object.
  m_PhotometricControlPoints = example.m_PhotometricControlPoints;
  m_PhotometricWeights = example.m_PhotometricWeights;

  m_PhotometricKernelType = example.m_PhotometricKernelType;
  m_PhotometricKernelWidth = example.m_PhotometricKernelWidth;

  m_PhotometricCPSpacing = example.m_PhotometricCPSpacing;
  m_NbPhotoCPs = example.m_NbPhotoCPs;
  m_NbVoxels = example.m_NbVoxels;

  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

  duplicator->SetInputImage(example.m_Image);
  duplicator->Update();
  m_Image = duplicator->GetOutput();

  duplicator->SetInputImage(example.m_DownSampledImage);
  duplicator->Update();
  m_DownSampledImage = duplicator->GetOutput();

  m_ImageSpatialGradient.resize(Dimension);
  for (unsigned int dim = 0; dim < Dimension; dim++) {
    duplicator->SetInputImage(example.m_ImageSpatialGradient[dim]);
    duplicator->Update();
    m_ImageSpatialGradient[dim] = duplicator->GetOutput(); }

  m_LIImage = example.m_LIImage;

  m_PermutationAxes = example.m_PermutationAxes;
  m_FlipAxes = example.m_FlipAxes;

  /// Deforms the image according to downSampledImageMap.
  if (example.m_DownSampledImage->GetLargestPossibleRegion().GetNumberOfPixels() != downSampledImageMap.rows())
    throw std::runtime_error("Number of voxels in downsampled image maps mismatch in ParametricImage constructor.");

  MatrixType Yup = GridFunctionsType::UpsampleImagePoints(m_Image, m_DownSampledImage, downSampledImageMap);
  InitializePhotometricKernelField(Yup);

  // Image.
  const VectorType I0 = m_PhotometricKernelField * m_PhotometricWeights;
  IteratorType it(m_Image, m_Image->GetLargestPossibleRegion());
  unsigned int k = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++k)
    it.Set(I0[k]);

  // Image spatial gradient.
  const ScalarType c(2 / pow(m_PhotometricKernelWidth, 2));
  ScalarType aux;
  for (unsigned int dim = 0 ; dim < Dimension ; ++dim)
  {
    it = IteratorType(m_ImageSpatialGradient[dim], m_ImageSpatialGradient[dim]->GetLargestPossibleRegion());
    k = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++k)
    {
      aux = 0;
      for (unsigned long l = 0 ; l < m_NbPhotoCPs ; ++l)
        aux += m_PhotometricKernelField(k, l) * m_PhotometricWeights(l)
            * (m_PhotometricControlPoints(l, dim) - Yup(k, dim));
      it.Set(c * aux);
    }
  }

  /// Finalization.
  this->SetParametricImageType();
  Update();
}


template <class ScalarType, unsigned int Dimension>
ParametricImage<ScalarType, Dimension>
::~ParametricImage()
{}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Encapsulation method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other public method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::Update()
{
  if (m_PhotometricWeights.size() == 0 || m_PhotometricControlPoints.rows() == 0)
    throw std::runtime_error("Photometric control points and weights should be set in ParametricImage.");

  if (this->IsModified())
  {
    this->UpdateBoundingBox();
    this->UnSetModified();
  }
}

template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::WriteObject(std::string str) const
{
  /// Writes on the native grid.
  VectorType I0 = GridFunctionsType::VectorizeImage(m_Image);
  m_LIImage->UpdateImageIntensity(I0);
  m_LIImage->WriteObject(str);

}

template <class ScalarType, unsigned int Dimension>
ScalarType
ParametricImage<ScalarType, Dimension>
::ComputeMatch(const std::shared_ptr<Superclass> target)
{
  if (!target->IsOfImageKind()) {
    std::cerr << "Exception : in ParametricImage::ComputeMatch(target), incompatible target deformable object."
              << std::endl; }

  Update();

  const std::shared_ptr<const LIImageType> LIItarget = std::static_pointer_cast<const LIImageType>(target);
  VectorType I0 = GridFunctionsType::VectorizeImage(m_Image);
  VectorType Ii = GridFunctionsType::VectorizeImage(LIItarget->GetImage());

  assert(I0.size() == Ii.size());

  VectorType D = I0 - Ii;
  return D.sum_of_squares();
}

// Computes nabla_{y(0)}A = 2 (I_0(y(0)) - I_i) nabla_{y(0)}I_0
template <class ScalarType, unsigned int Dimension>
MatrixType
ParametricImage<ScalarType, Dimension>
::ComputeMatchGradient(const std::shared_ptr<Superclass> target)
{
  if (!target->IsOfImageKind()) {
    std::cerr << "Exception : in ParametricImage::ComputeMatchGradient(target), incompatible target deformable object."
              << std::endl;
  }

  Update();

  const std::shared_ptr<const LIImageType> LIItarget = std::static_pointer_cast<const LIImageType>(target);
  VectorType I0 = GridFunctionsType::VectorizeImage(m_Image);
  VectorType Ii = GridFunctionsType::VectorizeImage(LIItarget->GetImage());
  VectorType D = I0 - Ii;

  MatrixType gradMatch(m_NbVoxels, Dimension);
  VectorType imageGrad_dim(m_NbVoxels);
  for (unsigned int dim = 0 ; dim < Dimension ; ++dim)
  {
    imageGrad_dim = GridFunctionsType::VectorizeImage(m_ImageSpatialGradient[dim]);
    gradMatch.set_column(dim, 2. * (D % imageGrad_dim));
  }

//    std::ostringstream oss;
//    oss << "0_ON_GradMatch.txt" << std::ends;
//    writeMatrixDLM<ScalarType>(oss.str().c_str(), gradMatch);

  return gradMatch;
}


template<class ScalarType, unsigned int Dimension>
MatrixType
ParametricImage<ScalarType, Dimension>
::SplatDifferenceImage(const std::shared_ptr<LIImageType> target, const MatrixType& downSampledImageMap)
{
  const MatrixType Yup = GridFunctionsType::UpsampleImagePoints(m_Image, m_DownSampledImage, downSampledImageMap);
  const ScalarType squaredKernelWidth = pow(m_PhotometricKernelWidth, 2);

  /// Instantiate the kernel object.
  KernelFactoryType* kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> ker = kfac->CreateKernelObject(m_PhotometricKernelType);
  ker->SetKernelWidth(m_PhotometricKernelWidth);

  /// Construct the photometric kernel field matrix.
  MatrixType photometricKernelField(m_NbVoxels, m_NbPhotoCPs, 0.0);
  for (unsigned int k = 0 ; k < m_NbVoxels ; ++k)
  {
    for (unsigned int p = 0 ; p < m_NbPhotoCPs ; ++p)
    {
      VectorType x = Yup.get_row(k);
      VectorType y = m_PhotometricControlPoints.get_row(p);
      photometricKernelField(k, p) = ker->EvaluateKernel(x, y);
    }
  }

  /// Compute the gradient wrt the photometric weights.
  const VectorType I0 = photometricKernelField * m_PhotometricWeights;
  const VectorType Ii = GridFunctionsType::VectorizeImage(target->GetImage());
  const VectorType D = I0 - Ii;

  return photometricKernelField.transpose() * D; // Equivalent of the splatting in the LinearInterpImage class.
}

template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::Initialize(std::shared_ptr<LIImageType> const LIImage)
{
  m_LIImage = LIImage;
  m_Image = LIImage->GetImage();
  m_DownSampledImage = LIImage->GetDownSampledImage();
  m_ImageSpatialGradient = LIImage->GetGradientImages();

  m_NbVoxels = m_Image->GetLargestPossibleRegion().GetNumberOfPixels();
  m_PermutationAxes = LIImage->GetPermutationAxes();
  m_FlipAxes = LIImage->GetFlipAxes();

  UpdateBoundingBox();
  InitializePhotometricControlPoints();
  InitializePhotometricWeights();

  const MatrixType Y0 = GridFunctionsType::ImageToPoints(m_Image);
  InitializePhotometricKernelField(Y0);
}


template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::UpdateImageIntensity(VectorType const& pw)
{
  assert(pw.size() == m_NbPhotoCPs);
  m_PhotometricWeights = pw;

  const MatrixType Y0 = GridFunctionsType::ImageToPoints(m_Image);
  InitializePhotometricKernelField(Y0);

  const VectorType I0 = m_PhotometricKernelField * m_PhotometricWeights;
  IteratorType it(m_Image, m_Image->GetLargestPossibleRegion());
  unsigned int k = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++k)
    it.Set(I0[k]);

  const ScalarType c(2 / pow(m_PhotometricKernelWidth, 2));
  ScalarType aux;
  for (unsigned int dim = 0 ; dim < Dimension ; ++dim)
  {
    it = IteratorType(m_ImageSpatialGradient[dim], m_ImageSpatialGradient[dim]->GetLargestPossibleRegion());
    k = 0;
    for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++k)
    {
      aux = 0;
      for (unsigned long l = 0 ; l < m_NbPhotoCPs ; ++l)
        aux += m_PhotometricKernelField(k, l) * m_PhotometricWeights(l)
            * (m_PhotometricControlPoints(l, dim) - Y0(k, dim));
      it.Set(c * aux);
    }
  }

  this->SetModified();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Protected method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::UpdateBoundingBox()
{
  MatrixType Y = GridFunctionsType::ImageToPoints(m_Image);

  for (unsigned int d = 0 ; d < Dimension ; d++)
  {
    Superclass::m_BoundingBox(d,0) = Y.get_column(d).min_value(); // - 0.5;
    Superclass::m_BoundingBox(d,1) = Y.get_column(d).max_value(); // + 0.5;
  }
}

template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::InitializePhotometricControlPoints()
{
  VectorType Xmin = Superclass::m_BoundingBox.get_column(0) - 1 * m_PhotometricKernelWidth;
  VectorType Xmax = Superclass::m_BoundingBox.get_column(1) + 1 * m_PhotometricKernelWidth;

  std::vector<VectorType> pointList;
  VectorType v(Dimension);
  switch (Dimension)
  {
    case 2:
    {
      ScalarType offsetX = 0.5*( Xmax[0] - Xmin[0] - m_PhotometricCPSpacing * floor( ( Xmax[0] - Xmin[0] ) / m_PhotometricCPSpacing ) );
      ScalarType offsetY = 0.5*( Xmax[1] - Xmin[1] - m_PhotometricCPSpacing * floor( ( Xmax[1] - Xmin[1] ) / m_PhotometricCPSpacing ) );

      for (ScalarType x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_PhotometricCPSpacing)
        for (ScalarType y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_PhotometricCPSpacing)
        {
          v[0] = x;
          v[1] = y;
          pointList.push_back(v);
        }
      break;
    }
    case 3:
    {
      ScalarType offsetX = 0.5*( Xmax[0] - Xmin[0] - m_PhotometricCPSpacing * floor( ( Xmax[0] - Xmin[0] ) / m_PhotometricCPSpacing ) );
      ScalarType offsetY = 0.5*( Xmax[1] - Xmin[1] - m_PhotometricCPSpacing * floor( ( Xmax[1] - Xmin[1] ) / m_PhotometricCPSpacing ) );
      ScalarType offsetZ = 0.5*( Xmax[2] - Xmin[2] - m_PhotometricCPSpacing * floor( ( Xmax[2] - Xmin[2] ) / m_PhotometricCPSpacing ) );

      for (ScalarType x = Xmin[0] + offsetX; x <= Xmax[0]; x += m_PhotometricCPSpacing)
        for (ScalarType y = Xmin[1] + offsetY; y <= Xmax[1]; y += m_PhotometricCPSpacing)
          for (ScalarType z = Xmin[2] + offsetZ; z <= Xmax[2]; z += m_PhotometricCPSpacing)
          {
            v[0] = x;
            v[1] = y;
            v[2] = z;
            pointList.push_back(v);
          }
      break;
    }
    default:
      throw std::runtime_error("ParametricImage::InitializePhotometricControlPoints is not "
                                   "implemented in Dimensions other than 2 and 3.");
  }

  m_NbPhotoCPs = pointList.size();

  m_PhotometricControlPoints.set_size(m_NbPhotoCPs, Dimension);
  for (unsigned int i = 0; i < m_NbPhotoCPs; i++)
    m_PhotometricControlPoints.set_row(i, pointList[i]);

  std::cout << "Set of " << m_NbPhotoCPs << " photometric control points defined." << std::endl;
}


template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::InitializePhotometricWeights()
{
  /// Instantiate the kernel object.
  KernelFactoryType* kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> ker = kfac->CreateKernelObject(m_PhotometricKernelType);
  ker->SetKernelWidth(m_PhotometricKernelWidth);

  /// Construct the Gram (symmetric) matrix G based on the photometric control points positions.
  MatrixType G(m_NbPhotoCPs, m_NbPhotoCPs);
  for (unsigned int i = 0 ; i < m_NbPhotoCPs ; ++i)
  {
    for (unsigned int j = 0 ; j < i ; ++j)
    {
      VectorType CPi = m_PhotometricControlPoints.get_row(i);
      VectorType CPj = m_PhotometricControlPoints.get_row(j);
      ScalarType k = ker->EvaluateKernel(CPi, CPj);

      G(i, j) = k;
      G(j, i) = k;
    }
    G(i, i) = 1;
  }

  /// Contruct the vector Ib of image values at the control points.
  VectorType Ib = GridFunctionsType::Interpolate(m_PhotometricControlPoints, m_LIImage->GetImage());

  /// Invert the linear system to find the weights.
  m_PhotometricWeights = solve(G, Ib);
}


template <class ScalarType, unsigned int Dimension>
void
ParametricImage<ScalarType, Dimension>
::InitializePhotometricKernelField(MatrixType const& grid)
{
  /// Instantiate the kernel object.
  KernelFactoryType* kfac = KernelFactoryType::Instantiate();
  std::shared_ptr<KernelType> ker = kfac->CreateKernelObject(m_PhotometricKernelType);
  ker->SetKernelWidth(m_PhotometricKernelWidth);

  /// Construct the photometric kernel field matrix.
  m_PhotometricKernelField.set_size(m_NbVoxels, m_NbPhotoCPs);
  for (unsigned int k = 0 ; k < m_NbVoxels ; ++k)
  {
    for (unsigned int p = 0 ; p < m_NbPhotoCPs ; ++p)
    {
      VectorType x = grid.get_row(k);
      VectorType y = m_PhotometricControlPoints.get_row(p);
      m_PhotometricKernelField(k, p) = ker->EvaluateKernel(x, y);
    }
  }
}

template class ParametricImage<double,2>;
template class ParametricImage<double,3>;
