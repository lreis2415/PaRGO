/***************************************************************************
* PeukerDouglas.cpp
*
* Project: GPRO_PD
* Purpose: Calculation of upward curved grid cells from DEM;Demonstration program for GPRO. 
*
* Author:  Fan Xingchen
* E-mail:  
****************************************************************************
* Copyright (c) 2022. Fan Xingchen
* 
****************************************************************************/


#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "PDoperator.h"
#include "utility.h"

using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: slope -elev <elevation grid file> -nbr <neighbor definition file> "
        "-out <output file>" << endl;

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
        } else { // Simple Usage
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
	demLayer.readNeighborhood(neighborfile); 
	demLayer.readFile(inputfilename,ROWWISE_DCMP);

	RasterLayer<double> ucgLayer("ucgLayer");//upward curved grid
	ucgLayer.copyLayerInfo(demLayer);
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	PDOperator pdOper;
	pdOper.demLayer(demLayer);
	pdOper.ucgLayer(ucgLayer);
	pdOper.Run();
	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;
	ucgLayer.writeFile(outputfilename);
	Application::END();
	return 0;
}
