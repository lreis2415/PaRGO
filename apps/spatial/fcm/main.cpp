/***************************************************************************
* main.cpp
*
* Project: PaRGO_0210
* Purpose: developing computation load balancing decomposition for GPRO. 
*			
* Author:  Zhan Lijun;Ai Beibei
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
#include "computeLayer.h"
#include "application.h"
#include "fcmOperator.h"
#include "communication.h"
#include "deComposition.h"

using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: fcm -elev <elevation grid file> " << endl
        << "-dataNbr <data layer neighbor file> " << endl
        << "-computeNbr <compute layer neighbor file>" << endl
        << "-out <output raster file>" << endl
        << "-clusterNum <number of clusters>" << endl
        << "-maxIter <max iteration number>" << endl
        << "-tolerance <tolerance in iteration>" << endl
        << "-weight <fuzzyness exponent>" << endl
        << "-dcmp <domain decomposition method>" << endl
        << "-writeLoad <path to write load file>" << endl
        << "-readLoad <path to read load file>" << endl << endl;

    cout << "'dcmp' available options:" << endl;
    cout << "\t space: (default) compute layer is decomposed equally by space, so it needs no evaluation." << endl;
    cout << "\t compute: compute layer is decomposed by computing load. To decide the load by running estimation before actual FCM algorithm"<< endl;

    cout << "(optional)'writeLoad' to create/rewrite time-cost load file to represent computing load. It"<< endl;
    cout << "\t either captures time cost if decompose set to 'space'" << endl;
    cout << "\t or estimates time cost if decompose set to 'compute'" << endl;

    cout << "(optional)'readLoad' is the path to load file to guide decomposition. Only needed when decompose set to 'compute'"<< endl;

    cout << " Or use the Simple Usage (full-text usage is recommended): slope <elevation grid file> "
        "<data layer neighbor file> "
        "<compute layer neighbor file> "
        "<output raster file> "
        "<number of clusters> "
        "<max iteration number> "
        "<tolerance in iteration> "
        "<fuzzyness exponent> "
        "[<domain decomposition method>] "
        "[<path to write load file>] "
        "[<path to read load file>] " << endl << endl;

    //cout << "Example.1. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif" << endl;
    //cout << "Example.2. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif -mtd SD" << endl;
    //cout << "Example.3. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif" << endl;
    //cout << "Example.4. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif TFD" << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
    /*  enum ProgramType{MPI_Type = 0,
                   MPI_OpenMP_Type,
                   CUDA_Type,
                   Serial_Type};*/
    /*  enum DomDcmpType{NON_DCMP = 0,
        ROWWISE_DCMP,
        COLWISE_DCMP,
        BLOCK_DCMP};*/
    /*  enum DomDcmpObj{SPACE_DIM = 0,
                   DATA_LOAD,
                   COMPT_LOAD};*/

    //...
    char* inputFileName;
    char* dataNeighbor;
    char* compuNeighbor;
    char* outputFileName;
    int clusterNum; //分类数目
    int maxIteration; //最大迭代次数
    double tolerance; //迭代阈值
    int weight; //加权指数
    bool decomposeBySapce; /// decomp by compute load if false
    char* writeLoadPath=nullptr;
    char* readLoadPath=nullptr;
    int i = 1;
    bool simpleusage = true;
    while (argc > i) {
        if (strcmp(argv[i], "-elev") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                inputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-elev'!");
            }
        }
        else if (strcmp(argv[i], "-dataNbr") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                dataNeighbor = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-dataNbr'!");
            }
        }
        else if (strcmp(argv[i], "-computeNbr") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                compuNeighbor = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-computeNbr'!");
            }
        }
        else if (strcmp(argv[i], "-out") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                outputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-out'!");
            }
        }
        else if (strcmp(argv[i], "-clusterNum") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                clusterNum = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-clusterNum'!");
            }
        }
        else if (strcmp(argv[i], "-maxIter") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                maxIteration = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-maxIter'!");
            }
        }
        else if (strcmp(argv[i], "-tolerance") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                tolerance = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-tolerance'!");
            }
        }
        else if (strcmp(argv[i], "-weight") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                weight = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-weight'!");
            }
        }
        else if (strcmp(argv[i], "-dcmp") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                decomposeBySapce = ifDecomposeBySpace(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-dcmp'!");
            }
        }
        else if (strcmp(argv[i], "-writeLoad") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                writeLoadPath = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-writeLoad'!");
            }
        }
        else if (strcmp(argv[i], "-readLoad") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                readLoadPath = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-readLoad'!");
            }
        }
        else {
            // Simple Usage
            if (!simpleusage) Usage("DO NOT mix the Full and Simple usages!");
            inputFileName = argv[1];
            dataNeighbor = argv[2];
            compuNeighbor = argv[3];
            outputFileName = argv[4];
            clusterNum = atoi(argv[5]);
            maxIteration = atoi(argv[6]);
            tolerance = atof(argv[7]);
            weight = atof(argv[8]);
            if (argc >= 10) {
                decomposeBySapce = ifDecomposeBySpace(argv[9]);
            }
            if(argc>=11){
                writeLoadPath = argv[10];
                readLoadPath = argv[10];
            }
            break;
        }
    }
    Application::START(MPI_Type, argc, argv); //init
    int myRank, process_nums;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_nums);

    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name,&name_len);
    if(myRank==0)
        cout<<"pargo-fcm. decompose by space: "<<decomposeBySapce<<". "<<process_nums<<" core(s)."<<endl;
    cout << "process " << myRank << " on " << processor_name << endl;

    double starttime;
    double endtime;

    //字符串解析输入文件名
    vector<const char *> vInputnames; //输入文件名，待解析
    vector<RasterLayer<double> *> vInputLayers;
    vector<string> tokens = SplitString(inputFileName,',');
    for(int i=0;i<tokens.size();i++) {
        vInputnames.push_back(tokens[i].c_str());
        RasterLayer<double>* pLayer = new RasterLayer<double>("none");
        vInputLayers.push_back(pLayer);
    }
    if (vInputnames.empty() || clusterNum == 0 || maxIteration == 0) {
        return 1;
    }
    RasterLayer<double> fcmLayer("fcmLayer"); //创建分类输出图层fcmLayer
    //预定义分类图层
    char** pDegLayerName = new char*[clusterNum];
    vector<RasterLayer<double> *> vDegreeLayer;
    for (int i = 0; i < clusterNum; i++) {
        pDegLayerName[i] = new char[50];
        sprintf(pDegLayerName[i], "degreeLayer%d.tif", i); //输出用名字
        RasterLayer<double>* pLayer = new RasterLayer<double>(pDegLayerName[i]);
        vDegreeLayer.push_back(pLayer);
    }

    const int comptGrain=10;
    if (decomposeBySapce) {
        //equal row dcmp based on region
        for (int i = 0; i < vInputnames.size(); i++) {
            vInputLayers[i]->readNeighborhood(dataNeighbor);
            vInputLayers[i]->readFile(vInputnames[i], ROWWISE_DCMP);
        }
        fcmLayer.copyLayerInfo(*vInputLayers[0]); //创建输出图层
        for (int i = 0; i < clusterNum; i++) {
            vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
        }

        FCMOperator fcmOper;
        fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, weight);
        fcmOper.inputLayer(vInputLayers);
        fcmOper.fcmLayer(fcmLayer);
        fcmOper.degLayer(vDegreeLayer);

        ComputeLayer<double> comptLayer("copmtLayer"); //暂时测试用，捕捉真实计算强度；以后改封装透明
        if(writeLoadPath) {
            comptLayer.init(vInputLayers[0],compuNeighbor);
            //comptLayer.copyLayerInfo(*vInputLayers[0]);
            //comptLayer.newMetaData(10);
            fcmOper.comptLayer(comptLayer);
        }
        starttime = MPI_Wtime();
        fcmOper.Run();
        if(writeLoadPath)
            comptLayer.writeFile(writeLoadPath); //测试用，写出捕捉到的计算时间
    }
    else{
        ////balanced row dcmp based on compute burden
        starttime = MPI_Wtime();
        vInputLayers[0]->readNeighborhood(dataNeighbor);
        CoordBR subWorkBR;
        ComputeLayer<double> comptLayer("computLayer");
        comptLayer.readNeighborhood(compuNeighbor);
        vInputLayers[0]->readGlobalFile(vInputnames[0]);
        comptLayer.addRasterLayer(*vInputLayers[0]);
        if(readLoadPath) {
            if(myRank==0) {
                comptLayer.readFile(readLoadPath);
                comptLayer.setComputGrain(comptGrain);
            }else {
                MPI_Barrier( MPI_COMM_WORLD );
            }
            comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
        }
        else {
            if (myRank == 0) {
                comptLayer.init(vInputLayers[0],compuNeighbor,comptGrain);
                Transformation<double> transOper(1, 15, &comptLayer); 
                transOper.run();
                comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
                if (writeLoadPath) {
                    comptLayer.writeComputeFile(writeLoadPath);               
                }
            } else {
                comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
                if(writeLoadPath) {
                    MPI_Barrier(MPI_COMM_WORLD);                    
                }
            }
        }
        cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows()
            << endl;
        endtime = MPI_Wtime();
        if (myRank == 0)
            cout << "dcmp time is " << endtime - starttime << endl;

        for (int i = 0; i < vInputnames.size(); i++) {
            vInputLayers[i]->readNeighborhood(dataNeighbor);
            vInputLayers[i]->readFile(vInputnames[i], subWorkBR, ROWWISE_DCMP);
        }
        fcmLayer.copyLayerInfo(*vInputLayers[0]); //创建输出图层
        for (int i = 0; i < clusterNum; i++) {
            vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
        }
        //执行计算
        FCMOperator fcmOper;
        fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, weight);
        fcmOper.inputLayer(vInputLayers);
        fcmOper.fcmLayer(fcmLayer);
        fcmOper.degLayer(vDegreeLayer);
        starttime = MPI_Wtime();
        fcmOper.Run();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();
    if (myRank == 0)
        cout << "compute time is " << endtime - starttime << endl<<endl;
    fcmLayer.writeFile(outputFileName);
    //for( size_t i = 0; i < vDegreeLayer.size(); ++i ){
    //	vDegreeLayer[i]->writeFile(pDegLayerName[i]);
    //}
    //cout << "write done." << endl;

    Application::END();
    return 0;
}
