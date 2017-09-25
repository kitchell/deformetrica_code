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

/// Class file.
#include "AbstractDeformations.h"

/// Core files.
#include "AdjointEquationsIntegrator.h"

/// Non-core files.
#include "itkImage.h"

/**
 *  \brief      Standard diffeomorphisms.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The Diffeos class inherited from AbstractDeformations represents the standard
 *              deformation which is usually employed in minimization problems such as registration
 *              or atlas construction. Deformations are encoded by control points and momentum vectors attached to them,
 *              which move in time according a to an Hamiltonian set of equations.\n \n
 */
template<class ScalarType, unsigned int Dimension>
class Diffeos : public AbstractDeformations<ScalarType, Dimension> {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Abstract deformation type.
  typedef AbstractDeformations<ScalarType, Dimension> Superclass;

  typedef DeformableMultiObject<ScalarType, Dimension> DeformableMultiObjectType;
  /// ITK image type.
  typedef itk::Image<ScalarType, Dimension> ImageType;
  /// ITK image pointer type.
  typedef typename ImageType::Pointer ImagePointerType;
  /// ITK image iterator type.
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIteratorType;
  /// Grid functions type.
  typedef GridFunctions<ScalarType, Dimension> GridFunctionsType;
  /// Kernel factory type.
  typedef KernelFactory<ScalarType, Dimension> KernelFactoryType;
  /// Kernel type.
  typedef typename KernelFactoryType::KernelBaseType KernelType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Contructor.
  Diffeos();

  /// Copy constructor.
  Diffeos(const Diffeos &other);

  /// Clones the object.
  std::shared_ptr<Diffeos> Clone() { return std::make_shared<Diffeos>(Diffeos(*this)); };

  /// Destructor.
  virtual ~Diffeos();

  // // Copy essential information, but not deformable multi-object, control points and momentas
  // virtual void CopyInformation(const Diffeos& other);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Return initial time.
  ScalarType GetT0() const { return m_T0; }
  /// Sets initial time to \e t0.
  void SetT0(ScalarType t0) {
    m_T0 = t0;
    this->SetModified();
  }

  /// Return final time.
  ScalarType GetTN() const { return m_TN; }
  /// Sets final time to \e tn.
  void SetTN(ScalarType tn) {
    m_TN = tn;
    this->SetModified();
  }

  /// Return the number of time points between \f$t_0\f$ and \f$t_n\f$.
  int GetNumberOfTimePoints() const { return m_NumberOfTimePoints; }
  /// Set the number of time points between \f$t_0\f$ and \f$t_n\f$.
  void SetNumberOfTimePoints(int n) {
    m_NumberOfTimePoints = n;
    this->SetModified();
  }

  void SetUseFastConvolutions() {
    m_UseFastConvolutions = true;
  }

  void UnsetUseFastConvolutions() {
    m_UseFastConvolutions = false;
  }

  ///	Return the type of the kernel.
  KernelEnumType GetKernelType() const { return m_KernelType; }
  /// Set the type of the kernel to \e kernelType.
  void SetKernelType(KernelEnumType kernelType) {
    m_KernelType = kernelType;
    this->SetModified();
  }

  ///	Return the size of the kernel.
  ScalarType GetKernelWidth() const { return m_KernelWidth; }
  /// Set the size of the kernel to \e kernelWidth.
  void SetKernelWidth(ScalarType kernelWidth) {
    m_KernelWidth = kernelWidth;
    this->SetModified();
  }

  /// Gets the initial control points.
  MatrixType GetStartPositions() const { return m_StartPositions; }
  /// Set the initial control points to \e X.
  void SetStartPositions(const MatrixType &X) {
    m_StartPositions = X;
    this->SetModified();
  }

  /// Set the initial momenta to \e A.
  void SetStartMomentas(const MatrixType &A) {
    m_StartMomentas = A;
    this->SetModified();
  }

  /// Get trajectories of control points
  MatrixListType GetTrajectoryPositions() const { return m_PositionsT; }
  /// Get values of momentum vectors in time
  MatrixListType GetTrajectoryMomentas() const { return m_MomentasT; }

  /// Returns the adjoint variable of CP positions (computed by solving adjoint equations).
  MatrixType GetAdjointPosAt0() const { return m_AdjointPosAt0; };
  /// Returns the adjoint variable of initial momenta (computed by solving adjoint equations).
  MatrixType GetAdjointMomAt0() const { return m_AdjointMomAt0; };
  /// Returns the adjoint variable of landmark points (computed by solving adjoint equations).
  MatrixType GetAdjointLandmarkPointsAt0() const { return m_AdjointLandmarkPointsAt0; };

  /// Return true if use improved Euler's method, false otherwise.
  bool ImprovedEuler() const { return m_UseImprovedEuler; }
  /// Set standard Euler's method.
  void UseStandardEuler() {
    m_UseImprovedEuler = false;
    this->SetModified();
  }
  /// Set improved Euler's method.
  void UseImprovedEuler() {
    m_UseImprovedEuler = true;
    this->SetModified();
  }

  /// Return the data domain.
  MatrixType GetDataDomain() const { return m_DataDomain; }
  /// Set the data domain to \e domain.
  void SetDataDomain(const MatrixType &domain) { m_DataDomain = domain; }

  /// Return true if any point is out of the bounding box, false otherwise.
  bool OutOfBox() const { return m_OutOfBox; }

  /// Return the padding factor.
  ScalarType GetPaddingFactor() const { return m_PaddingFactor; }
  /// Set the padding factor to \e paddingFactor.
  void SetPaddingFactor(ScalarType paddingFactor) { m_PaddingFactor = paddingFactor; }

  /// Return true if the true inverse flow is used, false otherwise.
  bool ComputeTrueInverseFlow() const { return m_ComputeTrueInverseFlow; }
  /// Set true inverse flow (see Diffeos::m_ComputeTrueInverseFlow for details).
  void SetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = true; }
  /// Set direct flow to compute inverse deformation (see Diffeos::m_ComputeTrueInverseFlow for details).
  void UnsetComputeTrueInverseFlow() { m_ComputeTrueInverseFlow = false; }

  /// Return true if the implicit Euler scheme is used, false otherwise.
  inline bool UseImplicitEuler() const { return m_UseImplicitEuler; }
  /// Sets the m_UseImplicitEuler flag to true.
  inline void SetUseImplicitEuler() { m_UseImplicitEuler = true; }
  /// Sets the m_UseImplicitEuler flag to false.
  inline void UnsetUseImplicitEuler() { m_UseImplicitEuler = false; }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Update();

  /// Computes the linear adjoint ODE of the flow equation (See AdjointEquationsIntegrator class for details).
  void IntegrateAdjointEquations(MatrixType &InitialConditionsLandmarkPoints,
                                 MatrixType &InitialConditionsImagePoints);
  void IntegrateAdjointEquations(MatrixListType &InitialConditionsLandmarkPoints,
                                 MatrixListType &InitialConditionsImagePoints,
                                 std::vector<unsigned int> jumpTimes);

  /**
   *	\brief		Checks if a set of points is out of box or not.
   *
   *	\details	This function enables to see if the coordinates of X at time t is out of box or not.
   *
   *	\param[in]	X	List of matrices containing coordinates at different time steps.
   *	\param[in]	t	Time index for X[t].
   *	\return		True if X[t] is out of box, false otherwise.
   */
  bool CheckBoundingBox(MatrixListType &X, int t);

// /*
//  *	\brief		Implements the linear adjoint ODE of the flow equation
//  *
//  *	\details	The flow equations writes: \f$\dot X(t) = G(X(t),S(t))\f$. This function computes \f$\dot\theta(t) = \partial_1 G(X(t),S(t))^T\theta(t) \f$ with \f$ \theta(1) \f$ the gradient of the fidelity term.
//  *
//  *	\param[in]  InitialConditions \f$ \theta(1) \f$	
//  *	\param[out]	PointsT	Trajectories of data and image points \f$ X(t) \f$
//  *	\param[out]	VectorsT Trajectories of auxiliary variable \f$ \theta(t) \f$
//  *	\return		True if X[t] is out of box, false otherwise.
//  */
// 	// void TransportAlongGeodesic(MatrixListType& InitialConditions, MatrixListType& VectorsT, MatrixListType& PointsT);

  virtual std::shared_ptr<DeformableMultiObjectType> GetDeformedObject() const {
    return this->GetDeformedObjectAt((this->GetNumberOfTimePoints() - 1));
  }

  /// Returns the deformed objects at time \e t.
  std::shared_ptr<DeformableMultiObjectType> GetDeformedObjectAt(unsigned int t) const;
  /// Returns the deformed control points at time \e t.
  MatrixType GetDeformedControlPointsAt(unsigned int t) const;

  virtual void WriteFlow(const std::vector<std::string> &name, const std::vector<std::string> &extension);

  /// Splat the residual (this deformed image - target image) defined on the final image map.
  MatrixType SplatResidualImage(const std::shared_ptr<DeformableMultiObjectType> target);

  /// Splats and sum the residuals (this deformed image - target images) defined on the final image map.
  MatrixType SplatResidualImages(const std::vector<std::shared_ptr<DeformableMultiObjectType>> targets,
                                 const std::vector<unsigned int> &timeIndices);

  /// Returns the parallel transported vector \e initialControlPoints, at times specified in the \e times vector.
  MatrixListType ParallelTransport(MatrixType const &initialMomenta,
                                   MatrixType const &initialControlPoints,
                                   ScalarType const &initialTime,
                                   std::vector<ScalarType> const &targetTimes,
                                   MatrixListType &velocities);
  /// Parallel transport, with some default values, option 1.
  MatrixListType ParallelTransport(MatrixType const &initialMomenta,
                                   std::vector<ScalarType> const &targetTimes) {
    MatrixListType aux;
    return ParallelTransport(initialMomenta, m_StartPositions, m_T0, targetTimes, aux);
  }
  /// Parallel transport, with some default values, option 2.
  MatrixListType ParallelTransport(MatrixType const &initialMomenta) {
    if (m_NumberOfTimePoints > 0) {
      std::vector<ScalarType> targetTimes(m_NumberOfTimePoints);
      const ScalarType h = (m_TN - m_T0) / (m_NumberOfTimePoints - 1);
      for (unsigned int k = 0; k < m_NumberOfTimePoints; ++k) { targetTimes[k] = m_T0 + k * h; }
      return ParallelTransport(initialMomenta, targetTimes);
    } else {
      return MatrixListType(0);
    }
  }

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Clones the object.
  virtual std::shared_ptr<Superclass> doClone() {
    return std::static_pointer_cast<Superclass>(std::make_shared<Diffeos>(Diffeos(*this)));
  };

  virtual void FlowImagePointsTrajectory();
  virtual void FlowLandmarkPointsTrajectory();

  /**
   *	\brief		Initializes the bounding box.
   *
   *	\details	This method initialize the Diffeos::m_BoundingBox attribute according to
   *				the data domain, the padding factor and the size of the kernel.
   */
  void InitBoundingBox();

  /// Solves the Hamiltonian system associated to the initial positions and momenta.
  void Shoot();

 private:

  /// Compute voxels trajectories using the direct flow integrated backward with speed flipped (i.e. \f$\phi_t\circ\phi_1^{-1}\f$).
  void IntegrateImagePointsBackward();
  /// Compute voxels trajectories using the flow of inverse deformations \f$\phi_t^{-1}\f$.
  void IntegrateImagePointsWithTrueInverseFlow();

  /// Explicit Euler step, main switch.
  MatrixType ComputeExplicitEulerStep(const unsigned int t, const MatrixType &v, const std::vector<unsigned int> &size);
  /// Explicit Euler step, order 41 in space.
  MatrixType ComputeExplicitEulerStep_1(const unsigned int t,
                                        const MatrixType &v,
                                        const std::vector<unsigned int> &size);
  /// Explicit Euler step, order 4 in space.
  MatrixType ComputeExplicitEulerStep_4(const unsigned int t,
                                        const MatrixType &v,
                                        const std::vector<unsigned int> &size);

  /// Implicit Euler step, main switch.
  void SolveImplicitEulerStep(const unsigned int t, const MatrixType &v, const std::vector<unsigned int> &size);
  /// Implicit Euler step, order 2 in space.
  void SolveImplicitEulerStep_2(const unsigned int t, const MatrixType &v, const std::vector<unsigned int> &size);

 protected :

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Initial time \f$ t_0 \f$.
  ScalarType m_T0;
  /// Final time \f$ t_n \f$.
  ScalarType m_TN;
  /// Number of time points between \f$ t_0 \f$ and \f$ t_n \f$.
  int m_NumberOfTimePoints;

  /// Matrix containing the initial control points (Size : N x Dimension).
  MatrixType m_StartPositions;
  /// Matrix containing the initial momenta (Size : N x Dimension).
  MatrixType m_StartMomentas;

  /// List containing the position of the control points at different time points.
  MatrixListType m_PositionsT;
  /// List containing the position of the momenta at different time points.
  MatrixListType m_MomentasT;

  /// Type of the kernel.
  KernelEnumType m_KernelType;
  /// Size of the kernel associated to the deformation.
  ScalarType m_KernelWidth;

  ///	Matrix containing the min (resp. the max) of the initial positions at first (resp. second) line.
  MatrixType m_DataDomain;
  ///	Box where any trajectory must not exit (Size : 2 x Dimension with the min (resp. the max) at first (resp. second) line).
  MatrixType m_BoundingBox;
  ///	Multiplier coefficient for the creation of the bounding box ("boundingBox = dataDomain +/- paddingFactor*kernelWidth").
  ScalarType m_PaddingFactor;
  ///	Boolean which prevents computations if any point (e.g. the coordinates of trajectory) is outside the bounding box.
  bool m_OutOfBox;

  /// Boolean which indicates if we use improved Euler's method or not.
  bool m_UseImprovedEuler;
  /// This parameter is used if there is an image.
  /// If set, true inverse flow will be used (i.e. \f$\phi_t^{-1}\f$).
  /// If not, direct flow will be integrated backward (speed flipped) to compute inverse deformation (i.e. \f$\phi_t\circ\phi_1^{-1}\f$).
  /// It moves particles backward from time t=1 back to time t=0 along the same trajectory as in the forward flow.
  /// In this case, results are different from inverse deformation if \f$t\f$ is different from \f$0\f$ and \f$1\f$.
  /// \warning For regression, always use true inverse flow.
  bool m_ComputeTrueInverseFlow;
  /// Use the more stable implicit Euler scheme for the inverse flow.
  bool m_UseImplicitEuler;
  /// Temporary flag, to be removed once the ComputeTrueInverseFlow integration scheme is fixed.
  bool m_RegressionFlag;
 public:
  void SetRegressionFlag() { m_RegressionFlag = true; }
 protected:

  /// Flow of voxel positions (with backward integration i.e. \f$\phi_t^{-1}\f$)
  MatrixListType m_MapsT;
  /// Flow of voxel positions (with true inverse flow i.e. \f$\phi_t\circ\phi_1^{-1}\f$).
  MatrixListType m_InverseMapsT;

  /// Velocity of landmark points
  MatrixListType m_LandmarkPointsVelocity;

  /// Trajectory of the whole vertices of Landmark type at different time steps.
  MatrixListType m_LandmarkPointsT;

  /// Adjoint variable of control points at time 0 (computed by solving adjoint equations)
  MatrixType m_AdjointPosAt0;

  /// Adjoint variable of momentum vectors at time 0  (computed by solving adjoint equations)
  MatrixType m_AdjointMomAt0;

  /// Adjoint variable of landmark points at time 0  (computed by solving adjoint equations)
  MatrixType m_AdjointLandmarkPointsAt0;

  /// Wether to use fast convolutions (experimental feature)
  bool m_UseFastConvolutions;

}; /* class Diffeos */

