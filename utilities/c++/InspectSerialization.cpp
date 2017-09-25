/***************************************************************************************
 *                                                                                      *
 *                                     Deformetrica                                     *
 *                                                                                      *
 *    Copyright Inria and the University of Utah. All rights reserved. This file is     *
 *    distributed under the terms of the Inria Non-Commercial License Agreement.        *
 *                                                                                      *
 *                                                                                      *
 ****************************************************************************************/


#include "DeformetricaConfig.h"
#include "SerializeDeformationState.h"

int main(int argc, char **argv) {

  if (argc != 2) {
    std::cerr << "Please, define one input filename for the ";
    return 0;
  }


  std::string filename(argv[1]);

  def::utils::DeformationState deformation_state;
  deformation_state.load(filename);

  std::cout << deformation_state;


    return 0;
}