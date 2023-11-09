/***************************************************************************
* dropanalysis.cpp
*
* Project: GPRO_Dropan
* Purpose: Calculate the T test for the drop of first order streams and higher order streams under given threshold; Demonstration program for GPRO. 
*
* Author:  Fan Xingchen
****************************************************************************
* Copyright (c) 2022. Fan Xingchen
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
#include "streamOperator.h"
#include "communication.h"
#include "transformation.h"
#include "computeLayer.h"
using namespace GPRO;

void Usage(const string& error_msg = "") {

    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout << " Usage: stream -net <strream net boolean grid file> -dir<d8 file> -nbr <neighbor definition file> -out <output streamOrder grid file>" << endl;
    cout << " Or use the Simple Usage: stream <strream net boolean grid file> <d8 file> <neighbor definition file> <output streamOrder grid file>" << endl << endl;
    cout << "Example.1. dropanalysis -net /path/to/net.tif -dir /path/to/d8.tif -nbr /path/to/moore.nbr -out /path/to/streamOrder.tif" << endl;
    cout << "Example.4. dropanalysis /path/to/net.tif /path/to/d8.tif /path/to/moore.nbr /path/to/streamOrder.tif" << endl;

    exit(1);
}

int main(int argc, char *argv[]) 
{

	/*!
	 * Parse input arguments.
	 * DO NOT start the application unless the required inputs are provided!
	 */
	if (argc < 5) {
        Usage("Too few arguments to run this program.");
	}
	// Input arguments
	//char* demfilename = nullptr;
	char* dirfilename = nullptr;
	//char* scafilename = nullptr;
	char* netfilename = nullptr;
	char* neighborfile = nullptr;
	char* outfile=nullptr;
	//double thresh;


    int i = 1;
    bool simpleusage = true;
	while (argc > i) {
		if (strcmp(argv[i], "-net") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                netfilename = argv[i];
                i++;
			} else {
				Usage("No argument followed '-net'!");
			}
		}else if (strcmp(argv[i], "-dir") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                dirfilename = argv[i];
                i++;
			} else {
				Usage("No argument followed '-dir'!");
			}
		}else if (strcmp(argv[i], "-nbr") == 0) {
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
                outfile = argv[i];
                i++;
			} else {
				Usage("No argument followed '-out'!");
			}
		}else { // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            //demfilename = argv[1];
            netfilename = argv[1];
            //scafilename = argv[3];
			dirfilename = argv[2];
			neighborfile = argv[3];
			outfile = argv[4];
			//thresh = atof(argv[5]);
			break;
        }
	}
	//if (!FileExists(demfilename)) {
	//Usage("The input DEM file not exists");
   // }
    if (!FileExists(neighborfile)) {
        Usage("neighbor file not exists");
    }
	
	Application::START(MPI_Type, argc, argv); //init
	int myRank, process_nums;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	RasterLayer<double> netLayer("netLayer"); //创建图层
	netLayer.readNeighborhood(neighborfile);
	netLayer.readFile(netfilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_dcmp

	RasterLayer<double> dirLayer("dirLayer"); //创建图层
	dirLayer.readNeighborhood(neighborfile);
	dirLayer.readFile(dirfilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_d
	
	//RasterLayer<double> scaLayer("scaLayer"); //创建图层
	//scaLayer.readNeighborhood(neighborfile);
	//scaLayer.readFile(scafilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_dcmp
	
	RasterLayer<double> outLayer("outLayer"); //创建图层
	outLayer.copyLayerInfo(netLayer);
	double opt;//t statistics
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	
	CoordBR subWorkBR;
	ComputeLayer<double> comptLayer("computLayer");
    comptLayer.addRasterLayerSerial(&netLayer);
	comptLayer.init(nullptr, 10);
    Transformation<double> transOper(0, 1, &comptLayer); 
    transOper.run();
    comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR); // Decompose the spatial computational domain.
	if (myRank == 0) cout << "dcmp time is " << MPI_Wtime() - starttime << endl;
    cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows() << endl;
	
	//cout<<"hh"<<endl;
	streamOperator DPOper;	
	DPOper.netLayer(netLayer);
	//DPOper.areaLayer(scaLayer);
	DPOper.dirLayer(dirLayer);
	//DPOper.demLayer(demLayer);
	DPOper.orderLayer(outLayer);
	DPOper.Run();
	
	

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;
	
	outLayer.writeFile(outfile);
	Application::END();
	//system("pause");
	return 0;
}