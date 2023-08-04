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
        "-slp <output slope file> [-mtd <algorithm>]" << endl;
    cout << "The available algorithm for slope are:" << endl;
    cout << "\t FD: (default) Third-order finite difference weighted by reciprocal of squared distance" << endl;
    cout << "\t FFD: Frame finite difference" << endl;
    cout << "\t MD: Maximum downslope" << endl;
    cout << "\t SD: Simple difference" << endl;
    cout << "\t SFD: Second-order finite difference" << endl;
    cout << "\t TFD: Third-order finite difference" << endl;
    cout << "\t TFDW: Third-order finite difference weighted by reciprocal of distance" << endl << endl;
    cout << " Or use the Simple Usage: slope <elevation grid file> <neighbor definition file> "
        "<output slope file> [<algorithm>]" << endl << endl;
    cout << "Example.1. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif" << endl;
    cout << "Example.2. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif -mtd SD" << endl;
    cout << "Example.3. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif" << endl;
    cout << "Example.4. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif TFD" << endl;

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
    SlopeAlgor calcalgor = FD;

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
        } else if (strcmp(argv[i], "-slp") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                outputfilename = argv[i];
                i++;
			} else {
                Usage("No argument followed '-slp'!");
			}
        } else if (strcmp(argv[i], "-mtd") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                calcalgor = GetSlopeAlgorithm(argv[i]);
                i++;
			} else {
                Usage("No argument followed '-mtd'!");
            }
        } else { // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputfilename = argv[1];
            neighborfile = argv[2];
            outputfilename = argv[3];
	        if (argc >= 5) {
                calcalgor = GetSlopeAlgorithm(argv[4]);
	        }
			break;
        }
	}

    if (!FileExists(inputfilename)) {
        Usage("The input DEM file not exists");
    }
    if (!FileExists(neighborfile)) {
        Usage("neighbor file not exists");
    }

	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type,
				   Serial_Type};*/
	Application::START(MPI_Type, argc, argv); //init

	RasterLayer<double> demLayer("demLayer");
	demLayer.readNeighborhood(neighborfile);
	demLayer.readFile(inputfilename,ROWWISE_DCMP);

	RasterLayer<double> slopeLayer("slopeLayer");
	slopeLayer.copyLayerInfo(demLayer);

	MPI_Barrier(MPI_COMM_WORLD);
	double starttime = MPI_Wtime();
	SlopeOperator slpOper;
	slpOper.demLayer(demLayer);
	slpOper.slopeLayer(slopeLayer);
	slpOper.calcAlgorithm(calcalgor);
	slpOper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	double endtime = MPI_Wtime();
	cout<<"run time is "<<endtime-starttime<<endl;

	slopeLayer.writeFile(outputfilename);
	
	Application::END();
	return 0;
}
