/***************************************************************************
* pTWI.cpp
*
* Project: GPRO_TWI
* Purpose: Calculation TWI from inputfile(SCA layer and slope layer);Demonstration program for GPRO. 
*
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2015. Ai Beibei
* 
****************************************************************************/


#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <sstream>
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "TWIOperator.h"
#include "mpi.h"
#include <omp.h>
#include "gdal_priv.h"

using namespace std;
using namespace GPRO;

int main(int argc, char *argv[]) 
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type};*/
	Application::START(MPI_Type, argc, argv); //init初始化，选择mpi版本

	char* inputfilename0;
	char* inputfilename1;
	char* neighborfile;
	char* outputfilename;

	if( argc!=5 ){
		cerr<<"please input right parameter.";
		return 0;
	}else{
		inputfilename0 = argv[1];	//slope layer
		inputfilename1 = argv[2];	//SCA layer
		neighborfile = argv[3]; 
		outputfilename = argv[4];
		//threadNUM = atoi(argv[5]);
	}

	RasterLayer<double> slopeLayer("slopeLayer");
	slopeLayer.readNeighborhood(neighborfile);
	slopeLayer.readFile(inputfilename0);

	RasterLayer<double> scaLayer("scaLayer");
	scaLayer.readNeighborhood(neighborfile);
	scaLayer.readFile(inputfilename1);

	RasterLayer<double> twiLayer("twiLayer");
	twiLayer.copyLayerInfo(slopeLayer);

	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	TWIOperator twiOper;
	
	twiOper.SCALayer(scaLayer);
	twiOper.slopeLayer(slopeLayer);
	twiOper.twiLayer(twiLayer);
	
	twiOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	cout<<"run time is "<<endtime-starttime<<endl;

	twiLayer.writeFile(outputfilename);

	Application::END();
	return 0;
}
