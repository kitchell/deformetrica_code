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

#include <fstream>
#include <sstream>
#include <string>

template<class ScalarType>
MatrixType readMatrixDLM(const char *fn) {

//	std::cout << "Reading " << fn << std::endl;

  std::ifstream infile(fn);

  if (!infile.is_open()) {
    std::cout << "Error opening file" << std::endl; // to re-route error stream to std::cout
    throw std::runtime_error("error opening file");
  }

  unsigned int numRows = 0;

  bool isfirst = true;
  std::string firstline;

  std::string line;
  while (!infile.eof()) {
    getline(infile, line);

    if (isfirst) {
      firstline = line;
      isfirst = false;
    }

    numRows++;
  }

  unsigned int numCols = 0;

  std::stringstream ss(firstline);
  while (!ss.eof()) {
    float f;
    ss >> f;
    numCols++;
  }

  //infile.clear();
  //infile.seekg(0);
  infile.close();
  infile.open(fn);

  --numRows;

//	std::cout << "Reading " << numRows << " x " << numCols << " matrix" << std::endl;
  MatrixType M(numRows, numCols, 0);

  ScalarType a = 0;
  for (unsigned int i = 0; i < numRows; i++) {
    for (unsigned int j = 0; j < numCols; j++) {
      infile >> a;
      M(i, j) = a;
    }
  }

  infile.close();

  return M;
}

template<class ScalarType>
std::vector<MatrixType> readMultipleMatrixDLM(const char *fn) {

//	std::cout << "Reading " << fn << std::endl;

  std::ifstream infile(fn);

  if (!infile.is_open()) {
    std::cout << "error opening file: " << std::string(fn) << std::endl; // to re-route error stream to std::cout
    throw std::runtime_error("error opening file");
  }

  std::string firstline;
  getline(infile, firstline);

  std::stringstream ss(firstline);
  unsigned int N, numRows, numCols;
  ss >> N;
  ss >> numRows;
  ss >> numCols;
  //std::cout << "Reading " << N << " times a " << numRows << " x " << numCols << " matrix." << std::endl;

  if (N >= std::numeric_limits<unsigned int>::max() - 1000
      || numRows >= std::numeric_limits<unsigned int>::max() - 1000
      || numCols >= std::numeric_limits<unsigned int>::max() - 1000) {
    return std::vector<MatrixType>(0);
  }

  getline(infile, firstline);

  std::vector<MatrixType> M(N);
  MatrixType Mn(numRows, numCols, 0);

  for (unsigned int n = 0; n < N; n++) {
    Mn.fill(0);
    ScalarType a = 0;

    for (unsigned int i = 0; i < numRows; i++) {
      for (unsigned int j = 0; j < numCols; j++) {
        infile >> a;
        Mn(i, j) = a;
      }
    }
    M[n] = Mn;
    getline(infile, firstline);
  }

  infile.close();

  return M;
}

template<class ScalarType>
void writeMatrixDLM(std::string fn, const MatrixType &M) {
  unsigned int numRows = M.rows();
  unsigned int numCols = M.columns();
  if (numRows == 0 || numCols == 0)
    return;

  std::ofstream outfile(fn);

  for (unsigned int i = 0; i < numRows; i++) {
    for (unsigned int j = 0; j < (numCols - 1); j++)
      outfile << M(i, j) << " ";
    outfile << M(i, numCols - 1) << std::endl;
  }

  outfile.close();
}

template<class ScalarType>
void writeMultipleMatrixDLM(std::string fn, MatrixListType const &M) {
  std::ofstream outfile(fn);

  unsigned int N = M.size();
  if (N == 0)
    outfile.close();

  unsigned int numRows = M[0].rows();
  unsigned int numCols = M[0].columns();

  outfile << N << " " << numRows << " " << numCols << std::endl << std::endl;

  for (unsigned int n = 0; n < N; n++) {
    MatrixType Mn = M[n];
    for (unsigned int i = 0; i < numRows; i++) {
      for (unsigned int j = 0; j < (numCols - 1); j++)
        outfile << Mn(i, j) << " ";
      outfile << Mn(i, numCols - 1) << std::endl;
    }
    outfile << std::endl;
  }

  outfile.close();
}

template <class ScalarType>
void printMatrix(std::string const name, MatrixType const &M){

  std::cout << "Printing matrix : " << name << std::endl;
  for (unsigned int i=0;i<M.rows();++i) {
    std::cout << "Row " << i << " ";
    for (unsigned int j = 0; j < M.cols(); ++j)
      std::cout <<  M(i, j) << " ; ";
    std::cout << std::endl;
  }

}

template MatrixType readMatrixDLM<double>(const char *fn);
template std::vector<MatrixType> readMultipleMatrixDLM<double>(const char *fn);
template void writeMatrixDLM<double>(std::string fn, const MatrixType &M);
template void writeMultipleMatrixDLM<double>(std::string fn,MatrixListType const &M);
template void printMatrix<double>(std::string const name, MatrixType const &M);