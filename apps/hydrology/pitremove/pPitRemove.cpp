/***************************************************************************
* pPitRemove.cpp
*
* Project: GPRO_PitRemove
* Purpose: Pit&flat removing from DEM;Demonstration program for GPRO. 
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
#include "pitRemoveOperator.h"
#include "mpi.h"
//#include<omp.h>


using namespace std;
using namespace GPRO;

int main(int argc, char *argv[]) 
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type};*/
	Application::START(MPI_Type, argc, argv); //parallel version

	char* demfilename;
	char* neighborfile;
	char* pitfilename;
//	int threadNUM;

	if( argc!=4 ){
		cerr<<"wrong input parameter";
		return 0;
	}else{
		demfilename = argv[1];
		neighborfile = argv[2]; 
		pitfilename = argv[3];
		//threadNUM = atoi(argv[4]);
	}
	//omp_set_num_threads(threadNUM);
	
	RasterLayer<double> demLayer("demLayer");
	demLayer.readNeighborhood(neighborfile);
	demLayer.readFile(demfilename);

	RasterLayer<double> wdemLayer("wdemLayer");
	wdemLayer.readNeighborhood(neighborfile);
	wdemLayer.copyLayerInfo(demLayer);
	
	double starttime1;
	double endtime1;
	MPI_Barrier(MPI_COMM_WORLD);	//等待每个线程都输出自己的行列数
	starttime1 = MPI_Wtime();

	PitRemoveOperator PitOper;
	PitOper.demLayer(demLayer);
	PitOper.wdemLayer(wdemLayer);
	PitOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime1 = MPI_Wtime();

	cout<<"run time is "<<endtime1-starttime1<<endl;
	wdemLayer.writeFile(pitfilename);
	
	Application::END();
	return 0;
}