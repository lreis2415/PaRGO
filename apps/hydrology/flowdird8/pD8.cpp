/***************************************************************************
* pD8.cpp
*
* Project: GPRO_D8
* Purpose: Calculation of D8 flow direction from DEM; Demonstration program for GPRO. 
*
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* 
****************************************************************************/


#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <sstream>
#include <omp.h>
#include "mpi.h"
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "d8Operator.h"
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

	char* inputfilename;
	char* neighborfile;
	char* outputfilename;
	int dirType;
	//int threadNUM;
	/*
	* "dirType" used to name the representation type of flow direction result.
	*  "256" means: If the maximum descent to several cells is the same, the direction of neighborhood will be accumulated. 
	* "8" means only value between "0-7" will output. For the same value of maximum descent of several neighbors, the last one will be chosed.
	*/
	if( argc!=5 ){
		cerr<<"please input right parameter."<<endl;
		return 0;
	}else{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilename = argv[3];
		dirType = atoi(argv[4]);
		//threadNUM = atoi(argv[5]);
	}
	//omp_set_num_threads(threadNUM);
	if( !( dirType==8 || dirType==256 ) ){
		cerr<<"please input right parameter.";
		return 0;
	}

	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(inputfilename);  //读取栅格数据

	RasterLayer<double> sflowLayer("slopeLayer");
	sflowLayer.copyLayerInfo(demLayer);
	
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	D8Operator SFlOper(dirType);	
	SFlOper.demLayer(demLayer);
	SFlOper.d8Layer(sflowLayer);
	SFlOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;

	sflowLayer.writeFile(outputfilename);
	
	Application::END();
	//system("pause");
	return 0;
}