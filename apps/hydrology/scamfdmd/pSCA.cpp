/***************************************************************************
* pSCA.cpp
*
* Project: GPRO_D8_SCA
* Purpose: Calculation of SCA based on MFDmd flow direction; Demonstration program for GPRO. 
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

#include "utility.h"
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "scaOperator.h"
#include "communication.h"


using namespace std;
using namespace GPRO;

int main(int argc, char *argv[]) 
{
    /*  enum ProgramType{MPI_Type = 0,
                   MPI_OpenMP_Type,
                   CUDA_Type,
                   Serial_Type};*/
    Application::START(MPI_Type, argc, argv); //init

    string inputfilename;
    char* neighborfile;
    char* outputfilename;
    //int threadNUM;
    if( argc!=4 ){
        cerr<<"please input right parameter.";
        return 0;
    }else{
        inputfilename = argv[1];
        neighborfile = argv[2]; 
        outputfilename = argv[3];
        //threadNUM = atoi(argv[5]);
    }
    //omp_set_num_threads(threadNUM);

    vector<string> weightfilenameVec = SplitString(inputfilename, ',');

    vector<RasterLayer<double> *> weightLayerVec;
    for(int i=0; i<8; i++)
    {
        RasterLayer<double> *tmpLayer=new RasterLayer<double>("unnamed");
        weightLayerVec.push_back(tmpLayer);
        weightLayerVec[i]->readNeighborhood(neighborfile);
        weightLayerVec[i]->readFile(weightfilenameVec[i].c_str());
    }
    
    RasterLayer<double> scaLayer("scaLayer");
    scaLayer.copyLayerInfo(*weightLayerVec[0]);
    
    double starttime;
    double endtime;
    MPI_Barrier(MPI_COMM_WORLD);
    starttime = MPI_Wtime();

    SCAOperator scaOper;    
    scaOper.mfdLayer(weightLayerVec);
    scaOper.scaLayer(scaLayer);
    scaOper.Run();

    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();
    cout<<"run time is "<<endtime-starttime<<endl;

    scaLayer.writeFile(outputfilename);
    
    Application::END();
    return 0;
}
