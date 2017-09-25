#include "vtkPoints.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkDelaunay3D.h"
#include "vtkSphereSource.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "vtkParametricEllipsoid.h"
//#include "vtkParametricSuperEllipsoid.h"
#include "vtkParametricFunctionSource.h"

#include "vtkDataSetSurfaceFilter.h"
#include "vtkTriangleFilter.h"
#include "vtkLoopSubdivisionFilter.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

#include <cmath>

typedef vnl_vector<double> VNLVectorType;
typedef vnl_matrix<double> VNLMatrixType;


int main(int argc, char** argv)
{

	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0] << " outfilename PointSet1 PointSet2 etc.." << std::endl;
		return -1;
	}

	int numShapes = argc - 2;
	char* outfilename = argv[1];

	std::vector< vtkSmartPointer<vtkPolyData> > shapes(numShapes);
	int TotalNumPts = 0;

	VNLVectorType center(3, 0.0);
	for (int i = 0; i < numShapes; i++)
	{
		vtkSmartPointer<vtkPolyDataReader> reader =
				vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(argv[i+2]);
		reader->Update();
		vtkSmartPointer<vtkPolyData> sh = reader->GetOutput();
		shapes[i] = sh;

		int numPts = sh->GetNumberOfPoints();
		TotalNumPts += numPts;
		for (int k = 0; k < numPts; k++)
		{
			double p[3];
			sh->GetPoint(k,p);
			for (int dim = 0; dim < 3; dim++)
				center(dim) += p[dim];
		}
	}
	center /= TotalNumPts;

	VNLMatrixType covM(3,3,0.0);
	for (int i = 0; i < numShapes; i++)
		for (int k = 0; k < shapes[i]->GetNumberOfPoints(); k++)
		{
			double p[3];
			shapes[i]->GetPoint(k,p);
			for (int r = 0; r < 3; r ++)
				for (int c = 0; c < 3; c++)
					covM(r,c) += (p[r] - center[r])*(p[c] - center[c]);
		}

	covM /= (TotalNumPts - 1);

	vnl_symmetric_eigensystem<double> eig(covM);

	VNLMatrixType R = eig.V;

	//vtkSmartPointer<vtkParametricSuperEllipsoid> pellipse =
	//  vtkSmartPointer<vtkParametricSuperEllipsoid>::New();
	vtkSmartPointer<vtkParametricEllipsoid> pellipse =
			vtkSmartPointer<vtkParametricEllipsoid>::New();
	pellipse->SetXRadius(1.5 * sqrt(eig.D(0, 0)));
	pellipse->SetYRadius(1.5 * sqrt(eig.D(1, 1)));
	pellipse->SetZRadius(1.5 * sqrt(eig.D(2, 2)));
	//pellipse->SetN1(0.8);
	//pellipse->SetN2(0.8);

	vtkSmartPointer<vtkParametricFunctionSource> psource =
			vtkSmartPointer<vtkParametricFunctionSource>::New();
	psource->SetParametricFunction(pellipse);
	psource->SetUResolution(30);
	psource->SetVResolution(30);
	psource->SetWResolution(30);
	psource->Update();

	vtkSmartPointer<vtkPolyData> templatePD = psource->GetOutput();
	for (unsigned int i = 0; i < templatePD->GetNumberOfPoints(); i++)
	{
		double x[3];
		templatePD->GetPoint(i, x);

		VNLMatrixType y(3, 1, 0.0);
		for (unsigned int dim = 0; dim < 3; dim++)
			y(dim, 0) = x[dim];

		y = R * y;
		for (unsigned int dim = 0; dim < 3; dim++)
			x[dim] = y(dim, 0) + center[dim];

		templatePD->GetPoints()->SetPoint(i, x);
	}

	std::cout << "Ellipsoid centered at " << center << std::endl;

	vtkSmartPointer<vtkPolyDataWriter> writer =
			vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(outfilename);
#if (VTK_MAJOR_VERSION <= 5)
	writer->SetInput(templatePD);
#else
	writer->SetInputData(templatePD);
#endif
	writer->Update();


}


