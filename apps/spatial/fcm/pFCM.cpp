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


    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName;
    char* dataNeighbor;
    char* compuNeighbor;
    char* outputFileName;
    int clusterNum;
    int maxIteration;
    double tolerance;
    int weight;
    int nodataLoad;
    int validLoad;
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
    
#ifdef _DEBUG
    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name,&name_len);
    if (myRank == 0) {
        cout<<"PaRGO-FCM. "<<process_nums<< " core(s) by "<<processor_name << endl;
        for (int i = 0; i < argc; ++i) {
            cout<<argv[i];
            if(argv[i][0]=='-') 
                cout<<" ";
            else
                cout<<endl;
        }
        cout<<endl;
    }
#endif

    double starttime;
    double endtime;

    vector<const char *> vInputnames;
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
    RasterLayer<double> fcmLayer("fcmLayer");
    char** pDegLayerName = new char*[clusterNum];
    vector<RasterLayer<double> *> vDegreeLayer;
    for (int i = 0; i < clusterNum; i++) {
        pDegLayerName[i] = new char[50];
        sprintf(pDegLayerName[i], "degreeLayer%d.tif", i);
        RasterLayer<double>* pLayer = new RasterLayer<double>(pDegLayerName[i]);
        vDegreeLayer.push_back(pLayer);
    }

    if (decomposeBySapce) {//equal-area strategy
        for (int i = 0; i < vInputnames.size(); i++) {
            vInputLayers[i]->newLocalNbrhood();
            vInputLayers[i]->readFile(vInputnames[i], ROWWISE_DCMP);
        }
        fcmLayer.copyLayerInfo(*vInputLayers[0]);
        for (int i = 0; i < clusterNum; i++) {
            vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
        }

        FCMOperator fcmOper;
        fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, weight);
        fcmOper.inputLayer(vInputLayers);
        fcmOper.fcmLayer(fcmLayer);
        fcmOper.degLayer(vDegreeLayer);

        ComputeLayer<double> comptLayer("copmtLayer");

        //if record execution time for preliminary experiment.
        //if true, the output file will be replaced by the time-recorded layer.
        if (writeLoadPath) {
            fcmOper._writePreExpLoad=true;
        }
        starttime = MPI_Wtime();
        fcmOper.Run();
#ifdef _DEBUG
        cout<< "rank" <<myRank<<" Membership degree compute time: "<<fcmOper.computeTimeExceptLastCell<<"s"<<endl;
        cout<< "rank" <<myRank<<" reduce time: "<<fcmOper.reduceTime<<"s"<<endl;
#endif
    }
    else{//the proposed load-balancing strategy

         //Fill the spatial computational domain. This is a serial procedure.
    	for (int i = 0; i < vInputnames.size(); i++) {
            vInputLayers[i]->readGlobalFileSerial(vInputnames[i]);
		}
        CoordBR subWorkBR;
        starttime = MPI_Wtime();
        if(readLoadPath) {//preliminary experiment mode, read the raster layer with recorded execution time.
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.init(vInputLayers,nullptr,granularity);
            comptLayer.readComputeLoadFile(readLoadPath);
            comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR ); // Decompose the spatial computational domain.
            if (myRank == 0)
                cout << "compute-read-dcmp time is " << MPI_Wtime() - starttime << endl;
        }
        else {//estimate function mode
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.init(vInputLayers,nullptr,granularity);
            Transformation<double> transOper(nodataLoad, validLoad, &comptLayer); 
            transOper.run();
            if (writeLoadPath) {
                comptLayer.writeComputeIntensityFile(writeLoadPath);               
            }
            comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR); // Decompose the spatial computational domain.
            if (myRank == 0)
                cout << "compute-write-dcmp time is " << MPI_Wtime() - starttime << endl;
        }
        cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows()
            << endl;

        // Decompose the data domain using the extent of spatial computational subdomains
        for (int i = 0; i < vInputnames.size(); i++) {
            vInputLayers[i]->readFile(vInputnames[i], subWorkBR, ROWWISE_DCMP);
        }
        fcmLayer.copyLayerInfo(*vInputLayers[0]);
        for (int i = 0; i < clusterNum; i++) {
            vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
        }
        FCMOperator fcmOper;
        fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, weight);
        fcmOper.inputLayer(vInputLayers);
        fcmOper.fcmLayer(fcmLayer);
        fcmOper.degLayer(vDegreeLayer);
        starttime = MPI_Wtime();
        fcmOper.Run();
#ifdef _DEBUG
        cout<< "rank" <<myRank<<" Membership degree compute time: "<<fcmOper.computeTimeExceptLastCell<<"s"<<endl;
        cout<< "rank" <<myRank<<" reduce time: "<<fcmOper.reduceTime<<"s"<<endl;
#endif

    }
    MPI_Barrier(MPI_COMM_WORLD);
    endtime = MPI_Wtime();
    if (myRank == 0)
        cout << "compute time is " << endtime - starttime << endl<<endl;

    fcmLayer.writeFile(outputFileName);
    cout << "write done." << endl;

    Application::END();
    return 0;
}
