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
#include "LEOperator.h"
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
	char* steplength;
	//char* outputfile2name;
	int threadNUM;
	if (argc < 6)
	{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilename = argv[3];
		steplength = argv[4];
	}
	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(inputfilename);  //读取栅格数据


    RasterLayer<double> demLayer2("demLayer2"); //创建图层
    demLayer2.readNeighborhood(neighborfile);  //读取分析窗口文件
    demLayer2.readFile(inputfilename);  //读取栅格数据

	RasterLayer<double> LELayer("LELayer");

	LELayer.copyLayerInfo(demLayer);

	double starttime;
	double endtime;

	MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();
    cout <<"program start:"<< endl;
	LEOperator LEOper;

	LEOper.demLayer(demLayer);
	LEOper.LELayer(LELayer);
	
    int steplengthint = atoi(steplength);
    for (int i = 0; i < LEOper.GetRowNum() - 2; i = i + steplengthint)
    {
        LELayer.copyLayerInfo(demLayer2);
        LEOper.SetCurrentScale(i);
        LEOper.Run();
        MPI_Barrier(MPI_COMM_WORLD);
        string fileNameStr = string(outputfilename);
        string fileNameAtCurrentScale = fileNameStr.insert(fileNameStr.find_last_of('.'), "_" + std::to_string((long long)i + 1)); // out.tif => out_1.tif
        MPI_Barrier(MPI_COMM_WORLD);
        LELayer.writeFile(fileNameAtCurrentScale.data());

        endtime = MPI_Wtime();
        cout << i << " run time is " << endtime - starttime << endl;
    }

	endtime = MPI_Wtime();
	cout<<"total run time is "<<endtime-starttime<<endl;
	
	Application::END();
	return 0;
}
