#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <sstream>
#include<omp.h>
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "SCAOperator.h"
#include "mpi.h"


using namespace std;
using namespace GPRO;

int main(int argc, char *argv[]) 
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type};*/
	Application::START(MPI_OpenMP_Type, argc, argv); //init

	char* inputfilename;
    char* inputfilename2;
	char* neighborfile;
	char* outputfilename;
	int kc_meth;
	float StepRatio;//step ratio
	int threadNUM;//thread number
	//if (argc < 9)
	{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilename = argv[3];
		kc_meth =atoi(argv[4]);
		StepRatio = atof(argv[5]);
		threadNUM = atoi(argv[6]);
        inputfilename2 = argv[7];
	}
	
	double starttime2;//for recording total runtime
	double endtime2;
	starttime2 = MPI_Wtime();
	omp_set_num_threads(threadNUM);
	RasterLayer<double> demLayer("demLayer");//input layer
	demLayer.readNeighborhood(neighborfile);
	demLayer.readFile(inputfilename);
	RasterLayer<double> tmpLayer("tmpLayer");//input layer
	tmpLayer.readNeighborhood(neighborfile);
	tmpLayer.readFile(inputfilename2);

	RasterLayer<double> SCALayer("SCALayer");//output layer
	SCALayer.copyLayerInfo(demLayer);

	double starttime1;//for recording computing time
	double endtime1;
	starttime1 = MPI_Wtime();
	SCAOperator SCAOper;
	
	SCAOper.demLayer(demLayer,kc_meth,StepRatio,threadNUM);
	SCAOper.SCALayer(SCALayer);
    SCAOper.tmpLayer(tmpLayer);
	
	SCAOper.Run();

	endtime1 = MPI_Wtime();

	SCALayer.writeFile(outputfilename);
	endtime2 = MPI_Wtime();

    cout<<"[DEBUG][TIMESPAN][IO]"<< (endtime2-starttime2)-(endtime1-starttime1)  << endl;
    cout<<"[DEBUG][TIMESPAN][COMPUTING]"<< endtime1-starttime1 << endl;
    cout<<"[DEBUG][TIMESPAN][TOTAL]"<< endtime2-starttime2 << endl;  

	Application::END();
	return 0;
}
