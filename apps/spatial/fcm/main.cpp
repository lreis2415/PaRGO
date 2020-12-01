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
    
    cout << " Usage: fcm -srcData <elevation grid file> " << endl
        << "-dataNbr <data layer neighbor file> " << endl
        << "-computeNbr <compute layer neighbor file>" << endl
        << "-out <output raster file>" << endl
        << "-clusterNum <number of clusters>" << endl
        << "-maxIter <max iteration number>" << endl
        << "-tolerance <tolerance in iteration>" << endl
        << "-weight <fuzzyness exponent>" << endl
        << "-dcmp <domain decomposition method>" << endl
        << "-writeLoad <path to write load file>" << endl
        << "-readLoad <path to read load file>" << endl
        << "-nodataLoad <load of NoData cells>" << endl
        << "-validLoad <load of nonempty cells>" << endl << endl;

    cout << "'dcmp' available options:" << endl;
    cout << "\t space: (default) compute layer is decomposed equally by space, so it needs no evaluation." << endl;
    cout << "\t compute: compute layer is decomposed by computing load. To decide the load by running estimation before actual FCM algorithm"<< endl;

    cout << "(optional)'writeLoad' to create/rewrite time-cost load file to represent computing load. It"<< endl;
    cout << "\t either captures time cost if decompose set to 'space'" << endl;
    cout << "\t or estimates time cost if decompose set to 'compute'" << endl;

    cout << "(optional)'readLoad' is the path to load file to guide decomposition. Only needed when decompose set to 'compute'"<< endl;

    cout << "(optional)'nodataLoad' is required when 'dcmp' set to compute and 'writeLoad' is switched on. It is the estimated load of NoData raster cells."<< endl;
    cout << "(optional)'validLoad' is required when 'dcmp' set to compute and 'writeLoad' is switched on. It is the estimated load of nonempty raster cells."<< endl;

    //cout << "Example.1. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif" << endl;
    //cout << "Example.2. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif -mtd SD" << endl;
    //cout << "Example.3. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif" << endl;
    //cout << "Example.4. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif TFD" << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName;
    char* dataNeighbor;
    char* compuNeighbor;
    char* outputFileName;
    int clusterNum; //分类数目
    int maxIteration; //最大迭代次数
    double tolerance; //迭代阈值
    int weight; //加权指数
    int nodataLoad; //
    int validLoad; //
    int granularity=10;
    bool decomposeBySapce; /// decomp by compute load if false
    char* writeLoadPath=nullptr;
    char* readLoadPath=nullptr;
    int i = 1;
    bool simpleusage = true;
    while (argc > i) {
        if (strcmp(argv[i], "-inputs") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                inputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-inputs'!");
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
        else if (strcmp(argv[i], "-nodataLoad") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                nodataLoad = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-nodataLoad'!");
            }
        }
        else if (strcmp(argv[i], "-validLoad") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                validLoad = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-validLoad'!");
            }
        }
        else if (strcmp(argv[i], "-granularity") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                granularity = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-granularity'!");
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
            // Simple Usage is not supported.
            if (!simpleusage) Usage();
        }
    }
    Application::START(MPI_Type, argc, argv); //init
    int myRank, process_nums;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_nums);

    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name,&name_len);
    if (myRank == 0) {
        cout<<"PaRGO-FCM. "<<process_nums<<" core(s)"<<endl;
        for (int i = 0; i < argc; ++i) {
            cout<<argv[i];
            if(argv[i][0]=='-') 
                cout<<" ";
            else
                cout<<endl;
        }
        cout<<endl;
    }
    //cout << "process " << myRank << " on " << processor_name << endl;

    double starttime;
    double endtime;

    //字符串解析输入文件名
    vector<const char *> vInputnames; //输入文件名，待解析
    vector<RasterLayer<double> *> vInputLayers;
    vector<string> tokens = SplitString(inputFileName,',');
    for(int i=0;i<tokens.size();i++) {
        vInputnames.push_back(tokens[i].c_str());
    	vector<string> fileNameFrags=SplitString(tokens[i].c_str(),'\\');
        RasterLayer<double>* pLayer = new RasterLayer<double>(fileNameFrags[fileNameFrags.size()-1]);
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

    if (decomposeBySapce) {
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
            comptLayer.initSerial(vInputLayers,compuNeighbor);
            fcmOper.comptLayer(comptLayer);
        }
        starttime = MPI_Wtime();
        fcmOper.Run();
        if(writeLoadPath)
            comptLayer.writeComputeIntensityFileSerial(writeLoadPath); //测试用，写出捕捉到的计算时间
        cout<< "rank" <<myRank<<" Membership degree compute time: "<<fcmOper.computeTimeExceptLastCell<<"s"<<endl;
    }
    else{
    	for (int i = 0; i < vInputnames.size(); i++) {
            vInputLayers[i]->readNeighborhoodSerial(dataNeighbor);
            vInputLayers[i]->readGlobalFileSerial(vInputnames[i]);
		}
        CoordBR subWorkBR;
        starttime = MPI_Wtime();
        if(readLoadPath) {
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.initSerial(vInputLayers,compuNeighbor,granularity);
            comptLayer.readComputeLoadFile(readLoadPath);
            comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
            if (myRank == 0)
                cout << "compute-read-dcmp time is " << MPI_Wtime() - starttime << endl;
        }
        else {
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.initSerial(vInputLayers,compuNeighbor,granularity);
            Transformation<double> transOper(nodataLoad, validLoad, &comptLayer); 
            transOper.run();
            if (writeLoadPath) {
                comptLayer.writeComputeIntensityFileSerial(writeLoadPath);               
            }
            comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
            if (myRank == 0)
                cout << "compute-write-dcmp time is " << MPI_Wtime() - starttime << endl;
        }
        cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows()
            << endl;
        if (myRank == 0)
            cout << "dcmp time is " << MPI_Wtime() - starttime << endl;

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
        cout<< "rank" <<myRank<<" Membership degree compute time: "<<fcmOper.computeTimeExceptLastCell<<"s"<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();
    if (myRank == 0)
        cout << "compute time is " << endtime - starttime << endl<<endl;

    fcmLayer.writeFile(outputFileName);
    // for( size_t i = 0; i < vDegreeLayer.size(); ++i ){
    // 	vDegreeLayer[i]->writeFile(pDegLayerName[i]);
    // }
    cout << "write done." << endl;

    Application::END();
    return 0;
}
