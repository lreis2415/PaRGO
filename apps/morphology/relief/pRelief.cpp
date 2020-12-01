/***************************************************************************
* pAspect.cpp
*
* Project: GPRO, v 1.0
* Purpose: Demonstration program for GPRO. 
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

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: relief -elev <elevation grid file> -nbr <neighbor definition file> -out <output relief file>" << endl;
    cout << " Or use the Simple Usage: relief <elevation grid file> <neighbor definition file> <output relief file>" << endl << endl;
    cout << "Example.1. relief -elev /path/to/elev.tif -nbr /path/to/moore.nbr -out /path/to/relief.tif" << endl;
    cout << "Example.2. relief /path/to/elev.tif /path/to/moore.nbr /path/to/relief.tif" << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
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

    int i = 1;
    bool simpleusage = true;
    while (argc > i) {
        if (strcmp(argv[i], "-elev") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                inputfilename = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-elev'!");
            }
        }
        else if (strcmp(argv[i], "-nbr") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                neighborfile = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-nbr'!");
            }
        }
        else if (strcmp(argv[i], "-out") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                outputfilename = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-out'!");
            }
        }
        else {
            // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputfilename = argv[1];
            neighborfile = argv[2];
            outputfilename = argv[3];
            break;
        }
    }

    Application::START(MPI_Type, argc, argv); //init
    RasterLayer<double> demLayer("demLayer"); //创建图层
    demLayer.readNeighborhood(neighborfile); //读取分析窗口文件
    demLayer.readFile(inputfilename); //读取栅格数据

    RasterLayer<double> reliefLayer("reliefLayer");

    reliefLayer.copyLayerInfo(demLayer);

    double starttime;
    double endtime;

    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();
    reliefOperator relOper;

    relOper.demLayer(demLayer);

    relOper.reliefLayer(reliefLayer);

    relOper.Run();

    MPI_Barrier(MPI_COMM_WORLD);

    endtime = MPI_Wtime();


    cout << "run time is " << endtime - starttime << endl;
    reliefLayer.writeFile(outputfilename);


    Application::END();
    return 0;
}
