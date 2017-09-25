/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef _SparseDiffeoWriter_h
#define _SparseDiffeoWriter_h

#include "Diffeos.h"

#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>

#include "LinearAlgebra.h"
using namespace def::algebra;
/**
 *  \brief      A text writer of SparseDiffeo.
 *
 *  \copyright  Inria and the University of Utah
 *  \version    Deformetrica 2.0
 *
 *  \details    The SparseDiffeoWriter class enables to save in a text file all the informations
 *              concerning the deformation namely the trajectory of the control point and the momentas,
 *              the size of the kernel, etc.
 */
template<class ScalarType, unsigned int Dimension>
class SparseDiffeoWriter {
 public:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // typedef :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  ///	Possible type of kernels.
  typedef enum { Exact, P3M } KernelEnumType;

  /// Deformation type.
  typedef Diffeos<ScalarType, Dimension> DiffeosType;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor(s) / Destructor :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  SparseDiffeoWriter();

  ~SparseDiffeoWriter();


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Encapsulation method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Sets the filename to \e fn.
  inline void SetFileName(const char *fn) { m_FileName = fn; }

  /// Sets the type of the kernel to \e type.
  void SetKernelType(const char *type);

  /// Sets the deformation to \e def.
  inline void SetDiffeos(DiffeosType *def) { m_Def = def; }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Other method(s) :
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Save in a file all the informations concerning the deformation.
  void Update();

 protected:

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Attribute(s)
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Entity representing the deformation.
  DiffeosType *m_Def;

  /// String containing the name of the file which be will saved thanks to the Update() method.
  const char *m_FileName;

  /// Type of the kernel.
  KernelEnumType m_KernelType;

};

#endif
