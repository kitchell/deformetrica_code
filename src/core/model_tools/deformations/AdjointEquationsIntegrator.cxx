
#include "AdjointEquationsIntegrator.h"
#include "SimpleTimer.h"

#include <cassert>
#include <itkDerivativeImageFilter.h>

#include "MatrixDLM.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor(s) / Destructor :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
AdjointEquationsIntegrator<ScalarType, Dimension>
::AdjointEquationsIntegrator() :
	m_IsLandmarkPoints(false), m_IsImagePoints(false), m_HasJumps(false), m_KernelObj1(NULL),
	m_KernelObj2(NULL), m_KernelObj3(NULL), m_KernelObj4(NULL), m_UseFastConvolutions(false)
{}



template <class ScalarType, unsigned int Dimension>
AdjointEquationsIntegrator<ScalarType, Dimension>
::~AdjointEquationsIntegrator()
{}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Other method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
AdjointEquationsIntegrator<ScalarType, Dimension>
::Update()
{
	m_XiPosT.resize(m_NumberOfTimePoints);
	m_XiMomT.resize(m_NumberOfTimePoints);
	m_ThetaT.resize(m_NumberOfTimePoints);
	m_EtaT.resize(m_NumberOfTimePoints);

	if (m_IsLandmarkPoints)
		this->IntegrateAdjointOfLandmarkPointsEquations();

	if (m_IsImagePoints)
	{
		if (m_ComputeTrueInverseFlow)
			this->IntegrateAdjointOfImagePointsBackward();
		else
			this->IntegrateAdjointOfImagePointsForward();

	}

	this->IntegrateAdjointOfDiffeoParametersEquations();
}


template <class ScalarType, unsigned int Dimension>
void
AdjointEquationsIntegrator<ScalarType, Dimension>
::IntegrateAdjointOfLandmarkPointsEquations()
{

	// Propagate theta backward. Initialization.
	m_ThetaT.resize(m_NumberOfTimePoints); // <--- Initialize each matrix at 0 ?

    ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);
    int subjIndex = m_JumpTimes.size() - 1;

	if (m_HasJumps) // Initial condition for regression.
	{
        if (m_JumpTimes[subjIndex] == m_NumberOfTimePoints - 1)
        {
            m_ThetaT[m_NumberOfTimePoints - 1] = m_ListInitialConditionsLandmarkPoints[subjIndex];
            --subjIndex; if (subjIndex < 0) subjIndex = 0;
        }
        else
        {
            m_ThetaT[m_NumberOfTimePoints - 1] = m_ListInitialConditionsLandmarkPoints[0];
            m_ThetaT[m_NumberOfTimePoints - 1].fill(0.0);
        }
	}
    else
	{
		m_ThetaT[m_NumberOfTimePoints - 1] = m_ListInitialConditionsLandmarkPoints[0];
	}

	typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(m_KernelType);
	kernelObj->SetKernelWidth(m_KernelWidth);

	for (long t = m_NumberOfTimePoints - 1 ; t > 0; --t)
	{

		kernelObj->SetSources(m_PosT[t]);
		kernelObj->SetWeights(m_MomT[t]);

		MatrixType dTheta = kernelObj->ConvolveGradient(m_LandmarkPointsT[t], m_ThetaT[t]);

		if (m_HasJumps)
		{
			if (t == m_JumpTimes[subjIndex])
			{
				dTheta += m_ListInitialConditionsLandmarkPoints[subjIndex]; // WARNING : problem here if simple euler ?
			}
		}

		m_ThetaT[t-1] = m_ThetaT[t] + dTheta * dt; // the plus is correct! dTheta should be negative.

		// Heun's method
		if (m_UseImprovedEuler)
		{

			kernelObj->SetSources(m_PosT[t-1]);
			kernelObj->SetWeights(m_MomT[t-1]);

			MatrixType dTheta2 = kernelObj->ConvolveGradient(m_LandmarkPointsT[t-1], m_ThetaT[t-1]);

			if (m_HasJumps)
			{
				if (t == m_JumpTimes[subjIndex])
				{
					dTheta2 += m_ListInitialConditionsLandmarkPoints[subjIndex];
					--subjIndex; if (subjIndex < 0) subjIndex = 0;
				}
			}

			m_ThetaT[t-1] = m_ThetaT[t] + (dTheta + dTheta2) * (dt * 0.5f);
		}
	}

}


// This computes \dot{\eta}(t) = -\partial_1 G(t) \eta(t)
template <class ScalarType, unsigned int Dimension>
void
AdjointEquationsIntegrator<ScalarType, Dimension>
::IntegrateAdjointOfImagePointsBackward()
{
	MatrixType& YT0 = m_ImagePointsT[0];

    // Propagate eta backward. Initialization.
    m_EtaT.resize(m_NumberOfTimePoints); // <--- Initialize each matrix at 0 ?

    ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);
    int subjIndex = m_JumpTimes.size() - 1;

    // The source term for integration. There is a + sign here contrary to in TransportAlongGeodesicForward.
    if (m_HasJumps) // Initial condition for regression.
    {
        if (m_JumpTimes[subjIndex] == m_NumberOfTimePoints - 1)
        {
            m_EtaT[m_NumberOfTimePoints - 1] = m_ListInitialConditionsImagePoints[subjIndex];
            --subjIndex; if (subjIndex < 0) subjIndex = 0;
        }
        else
        {
            m_EtaT[m_NumberOfTimePoints - 1] = m_ListInitialConditionsImagePoints[0];
            m_EtaT[m_NumberOfTimePoints - 1].fill(0.0);
        }
    }
    else
    {
        m_EtaT[m_NumberOfTimePoints - 1] = m_ListInitialConditionsImagePoints[0];
    }

	typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;

	typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
	typedef typename KernelFactoryType::KernelBaseType KernelType;
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(m_KernelType);
	kernelObj->SetKernelWidth(m_KernelWidth);

	// Integrate backwards in time
	for (long t = m_NumberOfTimePoints - 1 ; t > 0 ; --t)
	{
		kernelObj->SetSources(m_PosT[t]);
		kernelObj->SetWeights(m_MomT[t]);

		// The velocity is always v_t(y0)
		// MatrixType VtY0 = kernelObj->Convolve(m_DownSampledY1);
		MatrixType VtY0 = kernelObj->Convolve(YT0);

		// This gives us a Dim*Dim matrix at every pixel of the original image
		// std::vector<MatrixType> firstTerm = kernelObj->ConvolveGradient(m_DownSampledY1);

		std::vector<MatrixType> firstTerm;
		if (m_UseFastConvolutions)
			firstTerm = kernelObj->ConvolveGradientImageFast(YT0, m_DownSampledImage);
		else
			firstTerm = kernelObj->ConvolveGradient(YT0);



		// We need to have eta(t) in the form of an image, so we can compute the jacobian
		std::vector<ImageTypePointer> etaTImage;
		etaTImage.resize(Dimension);
		for (unsigned int d = 0; d < Dimension; d++)
		{
			// Create an image from each Dimension
			ImageTypePointer img = GridFunctionsType::VectorToImage(m_FullResolutionImage, m_EtaT[t].get_column(d));
			etaTImage[d] = img;
		}

		// Compute jacobian of eta(t)
		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
		std::vector<ImageTypePointer> gradImages;
		gradImages.resize(Dimension*Dimension);
		unsigned int indx = 0;
		// Compute the gradient in all directions
		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
		{
			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
			{
				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
				derivf->SetInput(etaTImage[dim1]);
				derivf->SetDirection(dim2);
				derivf->SetOrder(1);
				derivf->SetUseImageSpacingOn();
				derivf->Update();

				gradImages[indx] = derivf->GetOutput();
				indx++;
			}
		}

		// Now we have all the information we need to compute dEta, but we need to loop over all the pixels
		// and construct the 3x3 jacobian matrix and 3x1 v_t(y0)
		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;

		IteratorType it(m_FullResolutionImage, m_FullResolutionImage->GetLargestPossibleRegion());
		int matrixIndex = 0;
		int numVoxels = m_FullResolutionImage->GetLargestPossibleRegion().GetNumberOfPixels();
		MatrixType dEta(numVoxels, Dimension, 0);

		// Loop over the grid to construct jacobian matrices and compute dY
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			typedef typename ImageType::IndexType ImageIndexType;
			ImageIndexType imageIndex = it.GetIndex();
			MatrixType jacobian(Dimension, Dimension);

			unsigned int tempMatrixIndex = 0;
			// Build the (DimensionxDimension) jacobian matrix
			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
			{
				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
				{
					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
					tempMatrixIndex++;
				}
			}

			// The trace of the first term
			ScalarType traceValue = trace(firstTerm[matrixIndex]);

			// Get the (Dimension) vector corresponding to this points v_t(y0)
			VectorType VtY0k = VtY0.get_row(matrixIndex);

			dEta.set_row(matrixIndex, -traceValue*m_EtaT[t].get_row(matrixIndex) - jacobian*VtY0k);


			matrixIndex++;
		}


		if (m_HasJumps)
		{
			if (t == m_JumpTimes[subjIndex])
			{
				dEta += m_ListInitialConditionsImagePoints[subjIndex];
				--subjIndex;
				if (subjIndex < 0) subjIndex = 0;
			}
		}

		m_EtaT[t-1] = m_EtaT[t] + dEta * dt;
	}

	// We have computed Eta(t), the backwards propagator actually needs -(d_yp Y)^t Eta(t)
	for (long t = 0; t < m_NumberOfTimePoints ; ++t)
	{
		// We need to have Eta(t) in the form of an image, so we can compute the jacobian
		std::vector<ImageTypePointer> YTImage;
		YTImage.resize(Dimension);
		for (unsigned int d = 0; d < Dimension; d++)
		{
			// Create an image from each Dimension
			ImageTypePointer img = GridFunctionsType::VectorToImage(m_FullResolutionImage, m_ImagePointsT[t].get_column(d));
			YTImage[d] = img;
		}

		// Compute jacobian of eta(t)
		// Store the gradient images in a (Dimension*Dimension) vector dfx/dx, dfx/dy, dfx/dz, dfy/dx, dfy/dy ...
		std::vector<ImageTypePointer> gradImages;
		gradImages.resize(Dimension*Dimension);
		unsigned int indx = 0;
		// Compute the gradient in all directions
		for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
		{
			for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
			{
				typedef itk::DerivativeImageFilter<ImageType, ImageType> DerivativeFilterType;
				typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
				derivf->SetInput(YTImage[dim1]);
				derivf->SetDirection(dim2);
				derivf->SetOrder(1);
				derivf->SetUseImageSpacingOn();
				derivf->Update();

				gradImages[indx] = derivf->GetOutput();
				indx++;
			}
		}

		// Now we need to loop over all the pixels and construct the 3x3 jacobian matrix
		typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
		IteratorType it(m_FullResolutionImage, m_FullResolutionImage->GetLargestPossibleRegion());
		int matrixIndex = 0;

		// Loop over the grid to construct jacobian matrices and compute dY
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)
		{
			typedef typename ImageType::IndexType ImageIndexType;
			ImageIndexType imageIndex = it.GetIndex();
			MatrixType jacobian(Dimension, Dimension);

			unsigned int tempMatrixIndex = 0;
			// Build the (DimensionxDimension) jacobian matrix
			for (unsigned int dim1 = 0; dim1 < Dimension; dim1++)
			{
				for (unsigned int dim2 = 0; dim2 < Dimension; dim2++)
				{
					ImageTypePointer curGradImage = gradImages[tempMatrixIndex];
					jacobian(dim1, dim2) = curGradImage->GetPixel(imageIndex);
					tempMatrixIndex++;
				}
			}

			m_EtaT[t].set_row(matrixIndex, - jacobian.transpose() * m_EtaT[t].get_row(matrixIndex));
			matrixIndex++;
		}
	}

}



template <class ScalarType, unsigned int Dimension>
void
AdjointEquationsIntegrator<ScalarType, Dimension>
::IntegrateAdjointOfImagePointsForward()
{
//	/// For profiling.
//	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	/// Initialization of the size the eta maps.
	m_EtaT.resize(m_NumberOfTimePoints);

    /// Initialization of time-related variables.
	const ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);

	/// Initial source term.
    m_EtaT[0] = - m_ListInitialConditionsImagePoints[0];

    /// Kernel object instantiation.
	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();
	std::shared_ptr<KernelType> kernelObj = kFactory->CreateKernelObject(m_KernelType);
	kernelObj->SetKernelWidth(m_KernelWidth);

	for (unsigned int t = 0 ; t < m_NumberOfTimePoints - 1 ; ++t)
	{
		kernelObj->SetSources(m_PosT[t]);
		kernelObj->SetWeights(m_MomT[t]);

		// Grad mom splatted at CP locations, with convolutions evaluated at y(t-1)
		// This computes \eta_k^t alpha_p /nabla_1 K(y_k, c_p)
        MatrixType dEta;
        if (m_UseFastConvolutions)
          dEta = - kernelObj->ConvolveGradientImageFast(m_ImagePointsT[t], m_EtaT[t], m_DownSampledImage);
        else
          dEta = - kernelObj->ConvolveGradient(m_ImagePointsT[t], m_EtaT[t]);

        m_EtaT[t + 1] = m_EtaT[t] + dEta * dt;
		// Heun's method
		if (m_UseImprovedEuler)
		{
			kernelObj->SetSources(m_PosT[t + 1]);
			kernelObj->SetWeights(m_MomT[t + 1]);
            MatrixType dEta2;
            if (m_UseFastConvolutions)
			  dEta2 = kernelObj->ConvolveGradientImageFast(m_ImagePointsT[t + 1], m_EtaT[t + 1], m_DownSampledImage);
            else
              dEta2 = kernelObj->ConvolveGradient(m_ImagePointsT[t + 1], m_EtaT[t + 1]);
			m_EtaT[t + 1] = m_EtaT[t] + (dEta + dEta2) * (dt * 0.5f);
		}

	}

}




template <class ScalarType, unsigned int Dimension>
void
AdjointEquationsIntegrator<ScalarType, Dimension>
::IntegrateAdjointOfDiffeoParametersEquations()
{

	long numCP = m_PosT[0].rows();

	MatrixType zeroM(numCP, Dimension, 0);

	for (long t = 0; t < m_NumberOfTimePoints; t++)
	{
		m_XiPosT[t] = zeroM;
		m_XiMomT[t] = zeroM;
	}

	KernelFactoryType* kFactory = KernelFactoryType::Instantiate();

	m_KernelObj1 = kFactory->CreateKernelObject(m_KernelType);
	m_KernelObj1->SetKernelWidth(m_KernelWidth);
	m_KernelObj2 = kFactory->CreateKernelObject(m_KernelType);
	m_KernelObj2->SetKernelWidth(m_KernelWidth);
	m_KernelObj3 = kFactory->CreateKernelObject(m_KernelType);
	m_KernelObj3->SetKernelWidth(m_KernelWidth);

	ScalarType dt = (m_TN - m_T0) / (m_NumberOfTimePoints-1);

	for (long t = m_NumberOfTimePoints-1; t >= 1; t--)
	{
		MatrixType dPos;
		MatrixType dMom;
		this->ComputeUpdateAt(t, dPos, dMom);


		m_XiPosT[t-1] = m_XiPosT[t] + dPos*dt;
		m_XiMomT[t-1] = m_XiMomT[t] + dMom*dt;

		// Heun's method
		if (m_UseImprovedEuler)
		{
			MatrixType dPos2;
			MatrixType dMom2;
			this->ComputeUpdateAt(t-1, dPos2, dMom2);

			m_XiPosT[t-1] = m_XiPosT[t] + (dPos + dPos2) * (dt * 0.5f);
			m_XiMomT[t-1] = m_XiMomT[t] + (dMom + dMom2) * (dt * 0.5f);
		}

	} 
}








////////////////////////////////////////////////////////////////////////////////////////////////////
// Method(s) :
////////////////////////////////////////////////////////////////////////////////////////////////////

template <class ScalarType, unsigned int Dimension>
void
AdjointEquationsIntegrator<ScalarType, Dimension>
::ComputeUpdateAt(unsigned int s, MatrixType& dPos, MatrixType &dMom)
{

	long numCP = m_PosT[0].rows();

	std::shared_ptr<KernelType> momXiPosKernelObj = m_KernelObj1;
	std::shared_ptr<KernelType> etaKernelObj = m_KernelObj2;
	std::shared_ptr<KernelType> tmpKernelObj = m_KernelObj3;

	// Concatenate landmark and image points, as well as their adjoint variables. Save time in convolution.
	int nbOfLandmarkPoints = m_IsLandmarkPoints?m_LandmarkPointsT[0].rows():0;
	int nbOfImagePoints = m_IsImagePoints?m_ImagePointsT[0].rows():0;
	int nbTotalPoints = nbOfLandmarkPoints + nbOfImagePoints;

	MatrixType ConcatenatedPoints(nbTotalPoints, Dimension);
	MatrixType ConcatenatedVectors(nbTotalPoints, Dimension);

	for (int r = 0; r < nbOfLandmarkPoints; r++)
	{
		ConcatenatedPoints.set_row(r, m_LandmarkPointsT[s].get_row(r));
		ConcatenatedVectors.set_row(r, m_ThetaT[s].get_row(r));
	}
	for (int r = 0; r < nbOfImagePoints; r++)
	{
		// Be careful: if ComputeTrueInverseFlow, vectors m_EtaT[s] are always attached to the fixed points m_ImagePointsT[0]!
		ConcatenatedPoints.set_row(nbOfLandmarkPoints + r, m_ComputeTrueInverseFlow?m_ImagePointsT[0].get_row(r):m_ImagePointsT[s].get_row(r));
		ConcatenatedVectors.set_row(nbOfLandmarkPoints + r, m_EtaT[s].get_row(r));
	}

	etaKernelObj->SetSources(ConcatenatedPoints);
	etaKernelObj->SetWeights(ConcatenatedVectors);

	MatrixType dXi1 = etaKernelObj->ConvolveGradient(m_PosT[s], m_MomT[s]);

	MatrixType dXi2 = etaKernelObj->Convolve(m_PosT[s]);

	MatrixType AXiPos(numCP, Dimension*2, 0);
	AXiPos.set_columns(0, m_MomT[s]);
	AXiPos.set_columns(Dimension, m_XiPosT[s]);

	momXiPosKernelObj->SetSources(m_PosT[s]);
	momXiPosKernelObj->SetWeights(AXiPos);

	MatrixType kAXiPos = momXiPosKernelObj->Convolve(m_PosT[s]);

	MatrixListType gradAXiPos =	momXiPosKernelObj->ConvolveGradient(m_PosT[s]);

	MatrixType dXi3(numCP, Dimension, 0);
	for (unsigned int i = 0; i < numCP; i++)
	{
		MatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
		MatrixType gradXiPos_i = gradAXiPos[i].get_n_rows(Dimension, Dimension);
		dXi3.set_row(i,gradMom_i.transpose() * m_XiPosT[s].get_row(i) + gradXiPos_i.transpose() * m_MomT[s].get_row(i));
	}

	MatrixType dXi5 = kAXiPos.get_n_columns(Dimension, Dimension);

	MatrixType dXi6(numCP, Dimension, 0);

	for (unsigned int dim = 0; dim < Dimension; dim++)
	{
		MatrixType W(numCP, Dimension, 0.0);
		for (unsigned int i = 0; i < numCP; i++)
			W.set_row(i, m_XiMomT[s](i,dim) * m_MomT[s].get_row(i));

		tmpKernelObj->SetSources(m_PosT[s]);
		tmpKernelObj->SetWeights(W);
		MatrixType tmpgrad = tmpKernelObj->ConvolveGradient(m_PosT[s], dim);

		dXi6 += tmpgrad;
	}

	for (unsigned int i = 0; i < numCP; i++)
	{
		MatrixType gradMom_i = gradAXiPos[i].get_n_rows(0, Dimension);
		dXi6.set_row(i, dXi6.get_row(i) - gradMom_i * m_XiMomT[s].get_row(i));
	}

	tmpKernelObj->SetSources(m_PosT[s]);
	tmpKernelObj->SetWeights(m_MomT[s]);
	MatrixType dXi4 = tmpKernelObj->ConvolveSpecialHessian(m_XiMomT[s]);
	assert(dXi4.rows() == numCP);
	assert(dXi4.columns() == Dimension);


	dPos = dXi1 + dXi3 + dXi4;
	dMom = dXi2 + dXi5 + dXi6;

}




template class AdjointEquationsIntegrator<double,2>;
template class AdjointEquationsIntegrator<double,3>;


