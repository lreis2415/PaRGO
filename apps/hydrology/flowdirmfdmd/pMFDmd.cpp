/***************************************************************************
* pMFDmd.cpp
*
* Project: GPRO_MFDmd
* Purpose: Calculation of flow direction from DEM according to MFDmd (Qin2007); Demonstration program for GPRO. 
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
#include "mfdOperator.h"
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
	char* outputfilepath;
	double slpExp;
	//int threadNUM;
	/*
	* output the distribution weight matrix.
	*/
	if( argc!=5 ){
		cerr<<"please input right parameter.";
		return 0;
	}else{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilepath = argv[3];
		slpExp = atof(argv[3]);
		//threadNUM = atoi(argv[5]);
	}
	//omp_set_num_threads(threadNUM);
	if( !( slpExp>=0 || slpExp<=1 ) ){
		cerr<<"please input slope exponent between 0-1.";
		return 0;
	}

	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(inputfilename);  //读取栅格数据

	//creat output filenamme for weight matrix layers
	//string filepath = neighborfile;
	//size_t idx = filepath.find_last_of('/');
	//if( idx >= filepath.size() )
	//	idx = filepath.find_last_of('\\');
	//filepath.assign(filepath,0,idx+1);
	//const char* cFilepath = filepath.c_str();
	vector<RasterLayer<double> *> weightLayers;
	char** pOutputfilenames = new char*[8];
	for( int i=0; i<8; ++i ){
		char* tmp = new char [200];
		//sprintf(tmp, "%sweightLayer%d.tif", cFilepath,i);
		if( i<4 )
			sprintf(tmp, "%sweightLayer%d.tif", outputfilepath,i);
		else
			sprintf(tmp, "%sweightLayer%d.tif", outputfilepath,i+1);
		pOutputfilenames[i] = tmp;
		//RasterLayer<double> outLayer(pOutputfilenames[i]);
		RasterLayer<double> *outLayer = new RasterLayer<double> (pOutputfilenames[i]);
		weightLayers.push_back(outLayer);
		weightLayers[i]->copyLayerInfo(demLayer);
	}
	
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	MFDOperator SFlOper(slpExp);	
	SFlOper.demLayer(demLayer);
	SFlOper.weightLayers(weightLayers);
	SFlOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;

	for( int i=0; i<8; ++i ){
		weightLayers[i]->writeFile(pOutputfilenames[i]);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	//pOutputfilenames释放内存
	delete []pOutputfilenames;
	//weightLayers释放内存

	Application::END();
	return 0;
}