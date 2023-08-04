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
#include "dropanOperator.h"
#include "communication.h"

using namespace GPRO;

void Usage(const string& error_msg = "") {

    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout << " Usage: dropanalysis -elev <elevation grid file> -dir<d8 file> -sca<sca-d8 file> -ssa<weighted sca file> -nbr <neighbor definition file> -th <threshold>" << endl;
    cout << " Or use the Simple Usage: dropanalysis <elevation grid file> <d8 file> <sca-d8 file> <weighted sca file> <neighbor definition file> <threshold>" << endl << endl;
    cout << "Example.1. dropanalysis -elev /path/to/elev.tif -dir /path/to/d8.tif -sca /path/to/scad8.tif -ssa /path/to/w_sca.tif.tif -nbr /path/to/moore.nbr -th threshold" << endl;
    cout << "Example.4. dropanalysis /path/to/elev.tif /path/to/d8.tif /path/to/scad8.tif /path/to/w_sca.tif /path/to/moore.nbr threshold" << endl;

    exit(1);
}

int main(int argc, char *argv[]) 
{

	/*!
	 * Parse input arguments.
	 * DO NOT start the application unless the required inputs are provided!
	 */
	if (argc < 6) {
        Usage("Too few arguments to run this program.");
	}
	// Input arguments
	char* demfilename = nullptr;
	char* dirfilename = nullptr;
	char* scafilename = nullptr;
	char* ssafilename = nullptr;
	char* neighborfile = nullptr;
	double thresh;


    int i = 1;
    bool simpleusage = true;
	while (argc > i) {
		if (strcmp(argv[i], "-elev") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                demfilename = argv[i];
                i++;
            } else {
	            Usage("No argument followed '-elev'!");
            }
		} else if (strcmp(argv[i], "-dir") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                dirfilename = argv[i];
                i++;
			} else {
				Usage("No argument followed '-dir'!");
			}
		}else if (strcmp(argv[i], "-area") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                scafilename = argv[i];
                i++;
			} else {
				Usage("No argument followed '-area'!");
			}
		}else if (strcmp(argv[i], "-ssa") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                ssafilename = argv[i];
                i++;
			} else {
				Usage("No argument followed '-ssa'!");
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
		
        } else if (strcmp(argv[i], "-th") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                thresh = atof(argv[i]);
                i++;
			} else {
				Usage("No argument followed '-th'!");
			}
		}else { // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            demfilename = argv[1];
            dirfilename = argv[2];
            scafilename = argv[3];
			ssafilename = argv[4];
			neighborfile = argv[5];
			thresh = atof(argv[6]);
			break;
        }
	}
	if (!FileExists(demfilename)) {
        Usage("The input DEM file not exists");
    }
    if (!FileExists(neighborfile)) {
        Usage("neighbor file not exists");
    }
	
	Application::START(MPI_Type, argc, argv); //init

	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(demfilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_dcmp
	
	RasterLayer<double> ssaLayer("ssaLayer"); //创建图层
	ssaLayer.readNeighborhood(neighborfile);
	ssaLayer.readFile(ssafilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_dcmp

	RasterLayer<double> dirLayer("dirLayer"); //创建图层
	dirLayer.readNeighborhood(neighborfile);
	dirLayer.readFile(dirfilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_d
	
	RasterLayer<double> scaLayer("scaLayer"); //创建图层
	scaLayer.readNeighborhood(neighborfile);
	scaLayer.readFile(scafilename,ROWWISE_DCMP);  //读取栅格数据//add rowwise_dcmp
	
	
	double opt;//t statistics
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	
	
	//cout<<"thresh"<<thresh<<endl;
	DropanOperator DPOper(thresh,opt);	
		
	DPOper.ssaLayer(ssaLayer);
	DPOper.areaLayer(scaLayer);
	DPOper.dirLayer(dirLayer);
	DPOper.demLayer(demLayer);
	DPOper.Run();
	
	

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	//int myrank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//if( myrank==0 )
	cout<<"run time is "<<endtime-starttime<<endl;
	cout<<"t-statistics:"<<opt<<endl;
	Application::END();
	int result=2-fabs(opt);
	if(result>=0&&result<1)
	cout<<"optimal threshold:"<<thresh<<endl;
	//system("pause");
	return 0;
}