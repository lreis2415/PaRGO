/***************************************************************************
* pAspect.cpp
*
* Project: GPRO, v 1.0
* Purpose: Demonstration program for GPRO. 
*
* Usage:  
*
* Example: 
*
* Author:  Zhan Lijun
* E-mail:  zhanlj@lreis.ac.cn
****************************************************************************
* Copyright (c) 2013. Zhan Lijun
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
#include "reliefOperator.h"
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
	int threadNUM;
	if (argc < 6)
	{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilename = argv[3];
		//threadNUM = atoi(argv[4]);
	}
	//omp_set_num_threads(threadNUM);
	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(inputfilename);  //读取栅格数据

	RasterLayer<double> reliefLayer("reliefLayer");

	reliefLayer.copyLayerInfo(demLayer);
	

	double starttime;
	double endtime;

	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	reliefOperator relOper;
	
	//cout<<"SlopeOperator SlpOper;"<<endl;
	relOper.demLayer(demLayer);
	
	relOper.reliefLayer(reliefLayer);
	
	//cout<<"SlpOper.slopeLayer(slopeLayer);"<<endl;
	relOper.Run();
	
	MPI_Barrier(MPI_COMM_WORLD);

	endtime = MPI_Wtime();


	cout<<"run time is "<<endtime-starttime<<endl;
	reliefLayer.writeFile(outputfilename);
	
	
	Application::END();
	return 0;
}
