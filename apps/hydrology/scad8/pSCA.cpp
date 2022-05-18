/***************************************************************************
* pSCA.cpp
*
* Project: GPRO_D8_SCA
* Purpose: Calculation of SCA based on D8 flow direction; Demonstration program for GPRO. 
*
* Author:  Ai Beibei;Fan Xingchen
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


using namespace std;
using namespace GPRO;



int main(int argc, char* argv[]) {
    /*  enum ProgramType{MPI_Type = 0,
                   MPI_OpenMP_Type,
                   CUDA_Type,
                   Serial_Type};*/
	

    

    char* inputfilename;
    char* neighborfile;
    char* outputfilename;
	char* weightfile;//add weightfile
	bool usew;
    //int threadNUM;
    if (argc != 4&&argc != 5) {
        cerr << "please input right parameter.";
        return 0;
    }
    else if(argc==4){
        inputfilename = argv[1];
        neighborfile = argv[2];
        outputfilename = argv[3];
		weightfile="";
		usew=false;
        //threadNUM = atoi(argv[5]);
    }else{
		inputfilename = argv[1];
        neighborfile = argv[2];
		weightfile=argv[3];
        outputfilename = argv[4];
		usew=true;
		
	}
    //omp_set_num_threads(threadNUM);


	Application::START(MPI_Type, argc, argv); //init
    RasterLayer<double> d8Layer("d8Layer");
    d8Layer.readNeighborhood(neighborfile);
    d8Layer.readFile(inputfilename,ROWWISE_DCMP);//add rowwise_dcmp

    RasterLayer<double> scaLayer("scaLayer");
    scaLayer.copyLayerInfo(d8Layer);
	RasterLayer<double> weiLayer("weiLayer");
	if(usew){
		weiLayer.readNeighborhood(neighborfile);
		weiLayer.readFile(weightfile,ROWWISE_DCMP);
	}
	
		
	
    double starttime;
    double endtime;
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();

    SCAOperator scaOper;
    scaOper.d8Layer(d8Layer);
    scaOper.scaLayer(scaLayer);
	scaOper.usew=usew;
	if(usew)
		scaOper.weiLayer(weiLayer);

	
    scaOper.Run();

    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();
    cout << "run time is " << endtime - starttime << endl;

    scaLayer.writeFile(outputfilename);

    Application::END();
	//system("pause");
    return 10;
}
