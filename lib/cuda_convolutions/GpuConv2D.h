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

#ifndef _GpuConv2D_h
#define _GpuConv2D_h


////////////////////////////////////////////////////////////////////////////////////////////////////
// Conv2D :
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief   Computes the convolution of the kernel.
 *
 * \details This method computes \f$\gamma\f$ where :
 *          \f[
 *          \gamma_i = \sum_j K(x_i,y_j)\beta_j .
 *          \f]
 *
 * \param[in]   Ker     Evaluation kernel \f$K\f$.
 * \param[in]   x_h     Vector \f$x\f$ of size (nx*DIMPOINT).
 * \param[in]   y_h     Vector \f$y\f$ of size (ny*DIMPOINT).
 * \param[in]   beta_h  Vector \f$\beta\f$ of size (ny*DIMVECT).
 * \param[out]  gamma_h Vector \f$\gamma\f$ of size (nx*DIMVECT).
 * \param[in]   nx      Number of points of \f$x\f$.
 * \param[in]   ny      Number of points of \f$y\f$.
 */
template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuEvalConv2D(KER Ker, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);


// /// Performs GpuEvalConv2D with a Cauchy kernel.
// template < typename TYPE, int DIMPOINT, int DIMVECT >
// int CauchyGpuEvalConv2D(TYPE sigma, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);

/// Performs GpuEvalConv2D with a Gaussian kernel.
template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuEvalConv2D(TYPE sigma, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Grad1 Conv2D :
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief   Computes the gradient of the convolution.
 *
 * \details This method computes \f$\gamma\f$, the gradient w.r.t. the first variable of the convolution, where :
 *          \f[
 *          \gamma_i = 2 \sum_j (\alpha_i^T \beta_j) f'(r) \beta_j .
 *          \f]
 *
 * \param[in]   Ker     Evaluation kernel \f$K\f$.
 * \param[in]   x_h     Vector \f$x\f$ of size (nx*DIMPOINT).
 * \param[in]   alpha_h Vector \f$\alpha\f$ of size (nx*DIMVECT).
 * \param[in]   y_h     Vector \f$y\f$ of size (ny*DIMPOINT).
 * \param[in]   beta_h  Vector \f$\beta\f$ of size (ny*DIMVECT).
 * \param[out]  gamma_h Vector \f$\gamma\f$ of size (nx*DIMPOINT).
 * \param[in]   nx      Number of points of \f$x\f$.
 * \param[in]   ny      Number of points of \f$y\f$.
 */

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGrad1Conv2D(KER Ker, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);


// /// Performs GpuGrad1Conv2D with a Cauchy kernel.
// template < typename TYPE, int DIMPOINT, int DIMVECT >
// int CauchyGpuGrad1Conv2D(TYPE sigma, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);

/// Performs GpuGrad1Conv2D with a Gaussian kernel.
template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGrad1Conv2D(TYPE sigma, TYPE* alpha_h, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);


////////////////////////////////////////////////////////////////////////////////////////////////////
// Grad Conv2D :
////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 * \brief   Computes the convolution with the derivative of the Kernel (not to be confused with Grad1 which is the derivative of the convolution with the Kernel)
 *
 * \details This method computes \f$\gamma\f$, the convolution with the derivative (wrt to the \f$\ell\f$-th coordinate of \f$x\f$) of the Kernel, namely :
 *          \f[
 *          \gamma_i =  \sum_j \beta_j (f'(r))^{(\ell)}
 *          \f]
 *
 * \param[in]   Ker     Evaluation kernel \f$K\f$.
 * \param[in]   x_d     Vector \f$x\f$ of size (nx*DIMPOINT).
 * \param[in]   beta_d  Vector \f$\beta\f$ of size (nx*DIMVECT).
 * \param[in]   dim     Integer \f$\ell\f$ (0<dim<DIMPOINT).
 * \param[out]  gamma_d Vector \f$\gamma\f$ of size (nx*DIMVECT).
 * \param[in]   nx      Number of points of \f$x\f$.
 */

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGradConv2D(KER Ker,
        TYPE* x_h, TYPE* beta_h, unsigned int dim_h, TYPE* gamma_h,
         int nx);

///// Performs GpuGradConv2D with a Cauchy kernel.
//template < typename TYPE, int DIMPOINT, int DIMVECT >
//int CauchyGpuGradConv2D(TYPE sigma, TYPE* x_h, TYPE* beta_h, unsigned int dim_h, TYPE* gamma_h, int nx);

/// Performs GpuGradConv2D with a Gaussian kernel.
template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradConv2D(TYPE sigma, TYPE* x_h, TYPE* beta_h, unsigned int dim_h, TYPE* gamma_h, int nx);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Grad_varlin Conv2D :
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief   Computes the convolution with the derivative of the Kernel (not to be confused with Grad1 which is the derivative of the convolution with the Kernel)
 *
 * \details This method computes \f$\gamma\f$, the convolution with the derivative  of the Kernel, namely :
 *          \f[
 *          (\gamma_i)_{\ell} =  \sum_j \beta_j (f'(r))^{(\ell)}
 *          \f]
 *
 * \param[in]   Ker     Evaluation kernel \f$K\f$.
 * \param[in]   x_d     Vector \f$x\f$ of size (nx*DIMPOINT).
 * \param[in]   y_d     Vector \f$x\f$ of size (ny*DIMPOINT).
 * \param[in]   beta_d  Vector \f$\beta\f$ of size (nx*DIMVECT).
 * \param[out]  gamma_d Vector \f$\gamma\f$ of size (nx*(DIMVECT * DIMPOINT)).
 * \param[in]   nx      Number of points of \f$x\f$.
 * \param[in]   ny      Number of points of \f$y\f$.
 */

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGradConv_varlin_2D(KER Ker,
        TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h,
         int nx, int ny);

///// Performs GpuGradConv2D with a Cauchy kernel.
//template < typename TYPE, int DIMPOINT, int DIMVECT >
//int CauchyGpuGradConv2D(TYPE sigma, TYPE* alpha_h, TYPE* x_h, TYPE* beta_h, TYPE* gamma_h, int nx);

/// Performs GpuGradConv2D with a Gaussian kernel.
template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradConv_varlin_2D(TYPE sigma, TYPE* x_h, TYPE* y_h, TYPE* beta_h, TYPE* gamma_h, int nx, int ny);



////////////////////////////////////////////////////////////////////////////////////////////////////
// GradDiff Conv2D :
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief   Computes the hessian of the convolution .
 *
 * \details This method computes \f$\gamma\f$, the gradient of the differential
 *          along the vector\f$\eta\f$ of the convolution, where :
 *          \f[
 *          \gamma_i = 4 \sum_j (\beta_j^T \beta_i) f'(r) (\eta_i-\eta_j)
 *                   + 8 \sum_j(\beta_i^T \beta_j) f"(r) \left[(x_i-x_j)^T(\eta_i-\eta_j)\right](x_i-x_j)
 *          \f]
 *
 * \param[in]   Ker     Evaluation kernel \f$K\f$.
 * \param[in]   x_d     Vector \f$x\f$ of size (nx*DIMPOINT).
 * \param[in]   beta_d  Vector \f$\beta\f$ of size (nx*DIMVECT).
 * \param[in]   eta_d   Vector \f$\eta\f$ of size (nx*DIMPOINT).
 * \param[out]  gamma_d Vector \f$\gamma\f$ of size (nx*DIMPOINT).
 * \param[in]   nx      Number of points of \f$x\f$.
 */

template < typename TYPE, int DIMPOINT, int DIMVECT, class KER >
int GpuGradDiffConv2D(KER Ker,
        TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h,
         int nx);

template < typename TYPE, int DIMPOINT, int DIMVECT >
int GaussGpuGradDiffConv2D(TYPE sigma, TYPE* x_h, TYPE* beta_h, TYPE* eta_h, TYPE* gamma_h, int nx);


#endif /* _GpuConv2D_h */
