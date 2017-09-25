#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include <string>

#include "AnatomicalCoordinateSystem.h"

const unsigned int Dimension = 3;
typedef float ScalarType;

// typedef vnl_vector<ScalarType> VNLVectorType;
typedef vnl_matrix<ScalarType> VNLMatrixType;

typedef AnatomicalCoordinateSystem<ScalarType, Dimension> AnatomicalCoordinateSystemType;

using namespace std;

int main(int argc, char** argv)
{
	cout << "*************** Label 2 Matrix ***************" << endl;

	AnatomicalCoordinateSystemType ras_label2matrix("RAS");
	ras_label2matrix.PrintSelf();
	AnatomicalCoordinateSystemType sra_label2matrix("SRA");
	sra_label2matrix.PrintSelf();
	AnatomicalCoordinateSystemType lps_label2matrix("LPS");
	lps_label2matrix.PrintSelf();


	VNLMatrixType lps_matrix(Dimension, Dimension, 0.0);
	lps_matrix.set_identity();
	VNLMatrixType ras_matrix(Dimension, Dimension, 0.0);
	ras_matrix[0][0] = -1; ras_matrix[1][1] = -1; ras_matrix[2][2] = +1;
	VNLMatrixType sra_matrix(Dimension, Dimension, 0.0);
	sra_matrix[2][0] = +1; sra_matrix[0][1] = -1; sra_matrix[1][2] = -1;

	cout << "*************** Matrix 2 Label ***************" << endl;

	AnatomicalCoordinateSystemType ras_matrix2label(ras_matrix);
	ras_matrix2label.PrintSelf();
	AnatomicalCoordinateSystemType sra_matrix2label(sra_matrix);
	sra_matrix2label.PrintSelf();
	AnatomicalCoordinateSystemType lps_matrix2label(lps_matrix);
	lps_matrix2label.PrintSelf();

	if(ras_matrix2label.GetChangeOfBasisMatrix() != ras_label2matrix.GetChangeOfBasisMatrix())
		cerr << "Error with RAS !" << endl;
	if(sra_matrix2label.GetChangeOfBasisMatrix() != sra_label2matrix.GetChangeOfBasisMatrix())
		cerr << "Error with SRA !" << endl;
	if(lps_matrix2label.GetChangeOfBasisMatrix() != lps_label2matrix.GetChangeOfBasisMatrix())
		cerr << "Error with LPS !" << endl;

	return 0;
}

