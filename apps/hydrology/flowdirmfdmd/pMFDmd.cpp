/***************************************************************************
* pMFDmd.cpp
*
* Project: GPRO_MFDmd
* Purpose: Calculation of flow direction from DEM according to MFDmd (Qin et al., IJGIS, 2007).
*          A demonstration program for GPRO.
*
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* 
****************************************************************************/
#include <iostream>
#include <string>
#include <ctime>
#include "mpi.h"

#include "neighborhood.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "mfdOperator.h"

// using namespace std; // DO NOT use this statement in case of unexpected errors such as redefinitions of macros --LJ
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout <<
            " Usage: flowdirmfdmd -elev <elevation grid file> -nbr <neighbor definition file> -out <output flow direction file> -exp <slope exponent>"
            << endl;
    cout << " the default slope exponent is 1.1. For more details please see Qin 2007 or SimDTA manual" << endl << endl;
    cout <<
            " Or use the Simple Usage: flow <elevation grid file> <neighbor definition file> <output flow direction file> <slope exponent>"
            << endl << endl;
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
    double slpExp = -1.; // Should have a default value in case of no input --LJ
    char* slpExpStrErr = nullptr;

    int i = 1;
    bool simpleusage = true;
    while (argc > i) {
        if (strcmp(argv[i], "-elev") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                inputfilename = argv[i];
                i++;
            } else {
                Usage("No argument followed '-elev'!");
            }
        } else if (strcmp(argv[i], "-nbr") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                neighborfile = argv[i];
                i++;
            } else {
                Usage("No argument followed '-nbr'!");
            }
        } else if (strcmp(argv[i], "-out") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                outputfilename = argv[i];
                i++;
            } else {
                Usage("No argument followed '-out'!");
            }
        } else if (strcmp(argv[i], "-exp") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                // atof not report error, consider using strtod instead. --LJ
                // slpExp = atof(argv[i]);
                slpExp = strtod(argv[4], &slpExpStrErr);
                i++;
            } else {
                Usage("No argument followed '-exp'!");
            }
        } else {
            // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputfilename = argv[1];
            neighborfile = argv[2];
            outputfilename = argv[3];
            if (argc > 4) {
                // slpExp = atof(argv[4]);
                slpExp = strtod(argv[4], &slpExpStrErr);
            }
            break;
        }
    }
    if (!(slpExp >= 0 && slpExp <= 10)) {
        cerr << "please input slope exponent between 0-10.";
        return 0;
    }
	cout<<"slpexp:"<<slpExp<<endl;
    Application::START(MPI_Type, argc, argv); // init MPI environment

    RasterLayer<double> demLayer("demLayer"); // create raster layer (DEM)
    demLayer.readNeighborhood(neighborfile);  // read analysing neighbor window for raster layer
    demLayer.readFile(inputfilename, ROWWISE_DCMP);         // read data of layer (DEM)

    // creat output filenamme for weight matrix layers
    vector<RasterLayer<double>*> weightLayers;
    int flowdir_count = 8;
    char** pOutputfilenames = new char*[flowdir_count];
    for (int i = 0; i < flowdir_count; i++) {
        char* tmp = new char[200];
        // TODO: This is wrong.
        if (i < 4)
            sprintf(tmp, "%sweightLayer%d.tif", outputfilename, i);
        else
            sprintf(tmp, "%sweightLayer%d.tif", outputfilename, i + 1);
        pOutputfilenames[i] = tmp;
        RasterLayer<double>* outLayer = new RasterLayer<double>(pOutputfilenames[i]);
        weightLayers.push_back(outLayer);
        weightLayers[i]->copyLayerInfo(demLayer);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double starttime = MPI_Wtime();

    MFDOperator SFlOper(slpExp);
    SFlOper.demLayer(demLayer);
    SFlOper.weightLayers(weightLayers);
    SFlOper.Run();

    for (int i = 0; i < flowdir_count; i++) {
        weightLayers[i]->writeFile(pOutputfilenames[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double endtime = MPI_Wtime();

    int myrank = GetRank();
    // Should print the maximum time-consuming of all ranks. Data IO time should also be included. So, counting IO and computing separately. --LJ
    //if (myrank == 0)
    //    cout << "run time is " << endtime - starttime << endl;
    double t_all_rank = endtime - starttime;
    
    double t_all;  ///< Maximum time-consuming of parallel tasks in all ranks
    
    MPI_Reduce(&t_all_rank, &t_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        cout << "run time is " << t_all << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // For 2d-arrary, a for-loop should be used to release memory. --LJ
    // delete []pOutputfilenames;
    for (int i = 0; i < flowdir_count; i++) {
        delete[] pOutputfilenames[i];
        pOutputfilenames[i] = nullptr;
    }
    delete[] pOutputfilenames;
    pOutputfilenames = nullptr;

    Application::END();
    return 0;
}
