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

using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: pitremove -elev <elevation grid file> -nbr <neighbor definition file> -out <pit-removed output file>" << endl;
    cout << " Or use the Simple Usage: pitremove <elevation grid file> <neighbor definition file> <pit-removed output file>" << endl << endl;
    cout << "Example.1. pitremove -elev /path/to/elev.tif -nbr /path/to/moore.nbr -out /path/to/pitRemoved.tif" << endl;
    cout << "Example.2. pitremove /path/to/elev.tif /path/to/moore.nbr /path/to/pitRemoved.tif" << endl;

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


    Application::START(MPI_Type, argc, argv); //parallel version

    RasterLayer<double> demLayer("demLayer");
    demLayer.readNeighborhood(neighborfile);
    demLayer.readFile(inputfilename, ROWWISE_DCMP);

    RasterLayer<double> wdemLayer("wdemLayer");
    wdemLayer.readNeighborhood(neighborfile);
    wdemLayer.copyLayerInfo(demLayer);

    MPI_Barrier(MPI_COMM_WORLD);
    double starttime1 = 0;
    double endtime1 = 0;
    starttime1 = MPI_Wtime();

    PitRemoveOperator PitOper;
    PitOper.demLayer(demLayer);
    PitOper.wdemLayer(wdemLayer);
    PitOper.Run();

    MPI_Barrier(MPI_COMM_WORLD);

    endtime1 = MPI_Wtime();
    cout << "run time is " << endtime1 - starttime1 << endl;

    wdemLayer.writeFile(outputfilename);

    Application::END();
    return 0;
}
