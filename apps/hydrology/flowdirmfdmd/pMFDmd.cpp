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

void Usage(const string& error_msg = "") {

    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout << " Usage: flowdirmfdmd -elev <elevation grid file> -nbr <neighbor definition file> -out <output flow direction file> -exp <slope exponent>" << endl;
    cout << " the default slope exponent is 1.1. For more details please see Qin 2007 or SimDTA manual" << endl << endl;
    cout << " Or use the Simple Usage: flow <elevation grid file> <neighbor definition file> <output flow direction file> <slope exponent>" << endl << endl;
    cout << "Example.1. flowdirmfdmd -elev /path/to/elev.tif -nbr /path/to/moore.nbr -out /path/to/mfdmd.tif" << endl;
    cout << "Example.4. flowdirmfdmd /path/to/elev.tif /path/to/moore.nbr /path/to/mfdmd.tif" << endl;

    exit(1);
}

int main(int argc, char* argv[]) {

    /*!
     * Parse input arguments.
     * DO NOT start the application unless the required inputs are provided!
     */
    if (argc < 5) {
        Usage("Too few arguments to run this program.");
    }
    // Input arguments
    char* inputfilename = nullptr;
    char* neighborfile = nullptr;
    char* outputfilename = nullptr;
    double slpExp;

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
        else if (strcmp(argv[i], "-exp") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                slpExp = atof(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-exp'!");
            }
        }
        else {
            // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputfilename = argv[1];
            neighborfile = argv[2];
            outputfilename = argv[3];
            if (argc > 4) {
                slpExp = atof(argv[4]);
            }
            break;
        }
    }
    if (!(slpExp >= 0 || slpExp <= 1)) {
        cerr << "please input slope exponent between 0-1.";
        return 0;
    }

    Application::START(MPI_Type, argc, argv); //init

    RasterLayer<double> demLayer("demLayer"); //创建图层
    demLayer.readNeighborhood(neighborfile); //读取分析窗口文件
    demLayer.readFile(inputfilename); //读取栅格数据

    //creat output filenamme for weight matrix layers
    //string filepath = neighborfile;
    //size_t idx = filepath.find_last_of('/');
    //if( idx >= filepath.size() )
    //	idx = filepath.find_last_of('\\');
    //filepath.assign(filepath,0,idx+1);
    //const char* cFilepath = filepath.c_str();
    vector<RasterLayer<double>*> weightLayers;
    char** pOutputfilenames = new char*[8];
    for (int i = 0; i < 8; ++i) {
        char* tmp = new char [200];
        //sprintf(tmp, "%sweightLayer%d.tif", cFilepath,i);
        if (i < 4)
            sprintf(tmp, "%sweightLayer%d.tif", outputfilename, i);
        else
            sprintf(tmp, "%sweightLayer%d.tif", outputfilename, i + 1);
        pOutputfilenames[i] = tmp;
        //RasterLayer<double> outLayer(pOutputfilenames[i]);
        RasterLayer<double>* outLayer = new RasterLayer<double>(pOutputfilenames[i]);
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
    if (myrank == 0)
        cout << "run time is " << endtime - starttime << endl;

    for (int i = 0; i < 8; ++i) {
        weightLayers[i]->writeFile(pOutputfilenames[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //pOutputfilenames释放内存
    delete []pOutputfilenames;
    //weightLayers释放内存

    Application::END();
    return 0;
}
