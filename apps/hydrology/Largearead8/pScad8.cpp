/***************************************************************************
* pD8.cpp
*
* Project: GPRO_D8
* Purpose: Calculation of D8 flow direction from DEM; Demonstration program for GPRO. 
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
#include "scaOperator.h"
#include "communication.h"
#include "LargeareaProcess.h"
using namespace GPRO;
using namespace std;

void Usage(const string& error_msg = "") {

    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }
    cout << " Usage: flowdird8 -elev <elevation grid file> -nbr <neighbor definition file> -out <output flow direction file> -dirtype <direction type>" << endl;
    cout << " direction type is the representation type of flow direction result, should be one of the following:" << endl;
    cout << " \t 8: Only value between 0-7 will output. For the same value of maximum descent of several neighbors, the last one will be chosed." << endl;
    cout << " \t 256: If the maximum descent to several cells is the same, the direction of neighborhood will be accumulated." << endl;
    cout << " Or use the Simple Usage: flow <elevation grid file> <neighbor definition file> <output flow direction file>" << endl << endl;
    cout << "Example.1. flowdird8 -elev /path/to/elev.tif -nbr /path/to/moore.nbr -out /path/to/d8.tif" << endl;
    cout << "Example.4. flowdird8 /path/to/elev.tif /path/to/moore.nbr /path/to/d8.tif" << endl;

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
	char* inputwsfile = nullptr;
	char* inputfilename = nullptr;
	char* neighborfile = nullptr;
	char* outputfilename = nullptr;
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
		} else if (strcmp(argv[i], "-d") == 0) {
            simpleusage = false;
            i++;
			if (argc > i) {
                inputfilename = argv[i];
                i++;
            } else {
	            Usage("No argument followed '-d'!");
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
        } else { // Simple Usage
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
	//test(argc,argv);
	//init

	//if( !( dirType==8 || dirType==256 ) ){
		//cerr<<"please input right parameter.";
		//return 0;
	//}
	//cout<<"max:"<<argv[0];
	//RasterLayer<double> demLayer("demLayer"); //创建图层
	//读取分析窗口文件
	// //读取栅格数据//add rowwise_dcmp
	 
	//GDALAllRegister();
    //GDALDataset* poDatasetsrc = static_cast<GDALDataset*>(GDALOpen(inputfilename, GA_ReadOnly));

	
	//GDALRasterBand* poBandsrc = readFileInfo(poDatasetsrc, dcmpType);
	//poDatasetsrc->RasterIO(GF_Read, minICol, minIRow, nCols, nRows,
    //                    void *pData, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
    //                    _pMetaData->dataType, 0, 0);
	Application::START(MPI_Type, argc, argv);
	
	int pi,pj,qi,qj;
	LargeOperator LO;
	LO.srcfile=inputfilename;
	LO.outputfile=outputfilename;
	LO.wsfile=inputwsfile;
	LO.nbr=neighborfile;
	if(init){
		LO.initialfile();
		cout<<"create initial layer done"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,0);
		Application::END();
		return 0;
	}

	CoordBR sub=LO.getsubarea(pi,pj,qi,qj,_g,buf,id);
	
	/*
	//read watershed data and get the boundary coordinates of each sub-basin
	RasterLayer<double> iLayer("iLayer");
	iLayer.readFile(outputfilename,NON_DCMP);	
	RasterLayer<double> wsLayer("wsLayer");
	wsLayer.readFile(inputwsfile,NON_DCMP);
	Oper1.wslayer(wsLayer);
	Oper1.scaLayer(iLayer);
	Oper1.getarea(pi,pj,qi,qj,_g,buf,id);
	
	//create block of sub-basin
	CellCoord nw(pi,pj);
	CellCoord se(pi+qi-1,pj+qj-1);
	CoordBR sub(nw,se);
	*/
	//read sub-basin from whole flow-direction data
	RasterLayer<double> wdirLayer("wdirLayer");
	wdirLayer.readNeighborhood(neighborfile); 	 
	wdirLayer.readFile(inputfilename,sub,ROWWISE_DCMP);
	//read result of upstream watershed 
	RasterLayer<double> outLayer("outLayer");
	outLayer.readNeighborhood(neighborfile); 	 
	outLayer.readFile(outputfilename,sub,ROWWISE_DCMP);

	RasterLayer<double> scaLayer("scaLayer");
	scaLayer.copyLayerInfo(wdirLayer);
	
	double starttime;
	double endtime;
	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	
	SCAOperator d8Oper;
	d8Oper.outLayer(outLayer);
	d8Oper.d8Layer(wdirLayer);
	d8Oper.scaLayer(scaLayer);	
	d8Oper.Run();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	cout<<"run time is "<<endtime-starttime<<endl;

	//write out result of sub-basin
	string b="subwssca.tif";
	string a=to_string(_Longlong(id));
	string mid="data\\mid\\"+a+b;
	//string dem="dem\\"+a+"dem.tif";
	const char* p=mid.data();
	scaLayer.writeFile(p,sub);
	

	//write middle data into final data
	LO.pfile=p;
	LO.writefile(sub.nwCorner());
	/*
	//write middle data into final data
	RasterLayer<double> Layer("Layer");
	Layer.readNeighborhood(neighborfile); 	 
	Layer.readFile(p,NON_DCMP);
	Oper1.writesca(Layer,nw);
    iLayer.writeFile(outputfilename);
	*/
	Application::END();
	//system("pause");
	int time=int((endtime-starttime)*10000);
	return time;
	
}