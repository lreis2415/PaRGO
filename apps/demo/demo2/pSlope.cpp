/***************************************************************************
* pSlope.cpp
*
* Project: GPRO_Slope
* Purpose: Calculation of slope from DEM;Demonstration program for GPRO. 
* Review: Zhu Liangjun (07/18/2019 Reorganize the main function as a demo)
*
* Author:  Zhan Lijun; Ai Beibei; Wang Yijie
* E-mail:  zhanlj@lreis.ac.cn; aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2015-2019. Zhan Lijun;Ai Beibei
* 
****************************************************************************/


#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "slopeOperator.h"
#include "utility.h"

using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: slope -elev <elevation grid file> -nbr <neighbor definition file> "
        "-slp <output slope file>" << endl;

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
        else if (strcmp(argv[i], "-slp") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                outputfilename = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-slp'!");
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

    if (!FileExists(inputfilename)) {
        Usage("The input DEM file not exists");
    }
    if (!FileExists(neighborfile)) {
        Usage("neighbor file not exists");
    }

    Application::START(MPI_Type, argc, argv);

    RasterLayer<double> demLayer("demLayer");
    demLayer.readNeighborhood(neighborfile); // slope calculation on a raster cell needs its 3*3 neighborhood
    demLayer.readFile(inputfilename, ROWWISE_DCMP);

    RasterLayer<double> slopeLayer("slopeLayer");
    slopeLayer.copyLayerInfo(demLayer);

    SlopeOperator slpOper;
    slpOper.demLayer(demLayer);
    slpOper.slopeLayer(slopeLayer);
    slpOper.Run();

    slopeLayer.writeFile(outputfilename);

    Application::END();
    return 0;
}
