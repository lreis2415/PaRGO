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
#include "thOperator.h"
#include "communication.h"

using namespace GPRO;

void Usage(const string& error_msg = "") {

    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout << " Usage: flowdird8 -elev <elevation grid file> -nbr <neighbor definition file> -out <output flow direction file> -dirtype <direction type>" << endl;
    cout << " direction type is the representation type of flow direction result, should be one of the following:" << endl;
    cout << " \t 8: Only value between 0-7 will output. For the same value of maximum descent of several neighbors, the last one will be chosed." << endl;
    cout << " \t 256: If the maximum descent to several cells is the same, the direction of neighborhood will be accumulated." << endl;
    cout << " Or use the Simple Usage: flow <elevation grid file> <neighbor definition file> <output flow direction file>" << endl << endl;
    cout << "Example.1. flowdird8 -elev /path/to/elev.tif -nbr /path/to/moore.nbr -out /path/to/d8.tif" << endl;
    cout << "Example.4. flowdird8 /path/to/elev.tif /path/to/moore.nbr /path/to/d8.tif" << endl;

    exit(1);
}

int main(int argc, char *argv[]) 
{

	/*!
	 * Parse input arguments.
	 * DO NOT start the application unless the required inputs are provided!
	 */
	if (argc < 4) {
        Usage("Too few arguments to run this program.");
	}
	// Input arguments
	char* inputfilename = nullptr;
	char* neighborfile = nullptr;
	char* outputfilename = nullptr;
	double threshold=0.;

    int i = 1;
    bool simpleusage = true;
	while (argc > i) {
		if (strcmp(argv[i], "-sca") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                inputfilename = argv[i];
                i++;
            } else {
	            Usage("No argument followed '-sca'!");
            }
		}else if (strcmp(argv[i], "-net") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                outputfilename = argv[i];
                i++;
			} else {
                Usage("No argument followed '-net'!");
			}
        } else if (strcmp(argv[i], "-th") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                threshold = atof(argv[i]);
                i++;
			} else {
				Usage("No argument followed '-th'!");
			}
		}else { // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputfilename = argv[1];
            
            outputfilename = argv[2];
			threshold = atof(argv[3]);
			break;
        }
	}
	if (!FileExists(inputfilename)) {
        Usage("The input DEM file not exists");
    }
 
	
	Application::START(MPI_Type, argc, argv); //init


	RasterLayer<double> scaLayer("scaLayer"); //创建图层
	//scaLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	scaLayer.readFile(inputfilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_dcmp
	
	
	RasterLayer<double> netLayer("netLayer");
	netLayer.copyLayerInfo(scaLayer);
	
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();

	thOperator thOper(threshold);	
	thOper.scaLayer(scaLayer);
	thOper.netLayer(netLayer);
	thOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;
	netLayer.writeFile(outputfilename);
	Application::END();
	//system("pause");
	return 0;
}