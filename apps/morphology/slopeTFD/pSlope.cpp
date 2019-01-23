/***************************************************************************
* pSlope.cpp
*
* Project: GPRO_Slope
* Purpose: Calculation of slope from DEM;Demonstration program for GPRO. 
*
* Author:  Zhan Lijun;Ai Beibei
* E-mail:  zhanlj@lreis.ac.cn;aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2015. Zhan Lijun;Ai Beibei
* 
****************************************************************************/


#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <sstream>
#include<omp.h>
#include "mpi.h"
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "slopeOperator.h"
#include "communication.h"


using namespace std;
using namespace GPRO;

int main(int argc, char *argv[]) 
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type,
				   Serial_Type};*/
	Application::START(MPI_Type, argc, argv); //init

	//...
	char* inputfilename;
	char* neighborfile;
	char* outputfilename;
	//int threadNUM;

	if( argc != 4 ){
		cerr<<"please input right parameter";
	}else{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilename = argv[3];
		//threadNUM = atoi(argv[4]);
	}
	//omp_set_num_threads(threadNUM);
	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(inputfilename);  //读取栅格数据

	RasterLayer<double> slopeLayer("slopeLayer");
	slopeLayer.copyLayerInfo(demLayer);
	
	double starttime;
	double endtime;

	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	SlopeOperator slpOper;
	slpOper.demLayer(demLayer);
	slpOper.slopeLayer(slopeLayer);
	slpOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	cout<<"run time is "<<endtime-starttime<<endl;

	slopeLayer.writeFile(outputfilename);
	
	Application::END();
	return 0;
}
