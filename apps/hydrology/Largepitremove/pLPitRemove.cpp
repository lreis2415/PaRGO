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
#include "LpitRemoveOperator.h"
#include "LargeareaProcess.h"
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
	char* inputwsfile = nullptr;
	int _g=1;
	int buf=0;
	int id=0;
	bool init=false;

    int i = 1;
    bool simpleusage = true;
    while (argc > i) {
        if (strcmp(argv[i], "-w") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                inputwsfile = argv[i];
                i++;
            } else {
	            Usage("No argument followed '-ws'!");
            }
		} else if (strcmp(argv[i], "-elev") == 0) {
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
        }else if (strcmp(argv[i], "-g") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                _g = atoi(argv[i]);
                i++;
			} else {
                Usage("No argument followed '-g'!");
			}
        }else if (strcmp(argv[i], "-buf") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                buf = atoi(argv[i]);
                i++;
			} else {
                Usage("No argument followed '-buf'!");
			}
        }else if (strcmp(argv[i], "-id") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                id =atoi(argv[i]);
                i++;
			} else {
                Usage("No argument followed '-id'!");
			}
        }else {
            // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputwsfile = argv[1];
			inputfilename = argv[2];
            neighborfile = argv[3];
            outputfilename = argv[4];
			_g = atoi(argv[5]);
			buf = atoi(argv[6]);
			id = atoi(argv[7]);
            break;
        }
    }
	if (!FileExists(inputfilename)) {
        Usage("The input dir file not exists");
    }
    if (!FileExists(neighborfile)) {
        Usage("neighbor file not exists");
    }
	if (!FileExists(outputfilename)) {
        init=true;
    }

    Application::START(MPI_Type, argc, argv); //parallel version

	int pi,pj,qi,qj;
	LargeOperator LO;
	LO.srcfile=inputfilename;
	LO.outputfile=outputfilename;
	LO.wsfile=inputwsfile;
	LO.nbr=neighborfile;
	if(init){
		LO.initialfile();
		/*
		RasterLayer<double> dLayer("dLayer");
		dLayer.readNeighborhood(neighborfile);  
		dLayer.readFile(inputfilename, NON_DCMP);
		RasterLayer<double> initLayer("initLayer");
		initLayer.copyLayerInfo(dLayer);
		Oper1.init(initLayer);
		initLayer.writeFile(outputfilename);
		*/
		cout<<"create initial layer done"<<endl;
		//MPI_Finalize();
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,0);
		Application::END();
		return 0;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//read watershed data and get the boundary coordinates of each sub-basin
	CoordBR sub=LO.getsubarea(pi,pj,qi,qj,_g,buf,id);
	//cout<<"pi:"<<pi<<endl;
	/*
	RasterLayer<double> iLayer("iLayer");
	iLayer.readFile(outputfilename,NON_DCMP);	
	RasterLayer<double> wsLayer("wsLayer");
	wsLayer.readFile(inputwsfile,NON_DCMP);
	Oper1.wslayer(wsLayer);
	Oper1.wdemLayer(iLayer);
	Oper1.getarea(pi,pj,qi,qj,_g,buf,id);
	
	//create block of sub-basin
	CellCoord nw(pi,pj);
	CellCoord se(pi+qi-1,pj+qj-1);
	CoordBR sub(nw,se);
	*/
    RasterLayer<double> demLayer("demLayer");
    demLayer.readNeighborhood(neighborfile);
    demLayer.readFile(inputfilename, sub, ROWWISE_DCMP);
	//read result of upstream watershed 
	RasterLayer<double> outLayer("outLayer");
	outLayer.readNeighborhood(neighborfile); 	 
	outLayer.readFile(outputfilename,sub,ROWWISE_DCMP);

    RasterLayer<double> wdemLayer("wdemLayer");
    wdemLayer.readNeighborhood(neighborfile);
    wdemLayer.copyLayerInfo(demLayer);

    MPI_Barrier(MPI_COMM_WORLD);
    double starttime = 0;
    double endtime = 0;
    starttime = MPI_Wtime();

    PitRemoveOperator PitOper;
    PitOper.demLayer(demLayer);
    PitOper.wdemLayer(wdemLayer);
	PitOper.outLayer(outLayer);
    PitOper.Run();

    MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
    cout << "run time is " << endtime - starttime << endl;
	
	//write out result of sub-basin
	string b="subdemfel.tif";
	string a=to_string(_Longlong(id));
	string mid="data\\mid\\"+a+b;
	//string dem="dem\\"+a+"dem.tif";
	const char* p=mid.data();
	wdemLayer.writeFile(p,sub);
	MPI_Barrier(MPI_COMM_WORLD);

	//write middle data into final data

	LO.pfile=p;
	LO.writefile(sub.nwCorner());
	/*
	RasterLayer<double> Layer("Layer");
	Layer.readNeighborhood(neighborfile); 	 
	Layer.readFile(p,NON_DCMP);
	Oper1.writesca(Layer,nw);
	*/
	

    //iLayer.writeFile(outputfilename);
    Application::END();
	int time=int((endtime-starttime)*10000);
    return time;
}
