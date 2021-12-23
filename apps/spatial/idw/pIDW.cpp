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
#include <time.h>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computeLayer.h"
#include "application.h"
#include "idwOperator.h"
#include "transformation.h"
#include "idwTransformation.h"
using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: idw -sample <sample points file> " << endl
        << "-out <output raster file>" << endl
        << "-resolution <output resolution (length of one cell)(meter)>" << endl
        << "-fieldIndex <index of the field to be calculated>" << endl
        << "-idwExp <exponent of inverse distance>" << endl
        << "-searchPointNum <maximum of neighbor points to search>" << endl
        << "-searchRange <maximum of range to search>" << endl
        << "-blockSize <size of sample blocks>" << endl
        << "-dcmp <domain decomposition method>" << endl
        << "-writeLoad <path to write load file>" << endl
        << "-readLoad <path to read load file>" << endl << endl;

    cout << "'dcmp' available options:" << endl;
    cout << "\t space: (default) compute layer is decomposed equally by space, so it needs no evaluation." << endl;
    cout << "\t compute: compute layer is decomposed by computing load. To decide the load by running estimation before actual FCM algorithm" << endl;

    cout << "(optional)'writeLoad' to create/rewrite time-cost load file to represent computing load. It" << endl;
    cout << "\t either captures time cost if decompose set to 'space'" << endl;
    cout << "\t or estimates time cost if decompose set to 'compute'" << endl;

    cout << "(optional)'readLoad' is the path to load file to guide decomposition. Only needed when decompose set to 'compute'" << endl;


    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName=nullptr;
    char* maskFileName=nullptr;
    char* outputFileName=nullptr;
    float cellSize;
    int fldIdx;
    int idw_nbrPoints;
    int idw_power;
    double idw_buffer;
    double blockSize; //unit: meter
    int granularity = 10; //resolution of the computational domain = granularity * resolution of the data domain (10 by default)
    bool decomposeBySapce; /// decomp by compute load if false
    char* writeLoadPath = nullptr;
    char* readLoadPath = nullptr;

    int i = 1;
    bool simpleusage = true;
    while (argc > i) {
        if (strcmp(argv[i], "-sample") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                inputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-sample'!");
            }
        }
        if (strcmp(argv[i], "-mask") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                maskFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-mask'!");
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
        else if (strcmp(argv[i], "-resolution") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                cellSize = atof(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-resolution'!");
            }
        }
        else if (strcmp(argv[i], "-fieldIndex") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                fldIdx = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-fieldIndex'!");
            }
        }
        else if (strcmp(argv[i], "-idwExp") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                idw_power = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-idwExp'!");
            }
        }
        else if (strcmp(argv[i], "-searchPointNum") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                idw_nbrPoints = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-searchPointNum'!");
            }
        }
        else if (strcmp(argv[i], "-searchRange") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                idw_buffer = atof(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-searchRange'!");
            }
        }
        else if (strcmp(argv[i], "-blockSize") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                blockSize = atof(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-blockSize'!");
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
            // Simple Usage
            if (!simpleusage) Usage();
        }
    }

    Application::START(MPI_Type, argc, argv);
    int myRank, process_nums;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
    
#ifdef _DEBUG
    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);
    if (myRank == 0) {
        cout << "PaRGO-IDW. " << process_nums << " core(s) by "<<processor_name << endl;
        for (int i = 0; i < argc; ++i) {
            cout << argv[i];
            if (argv[i][0] == '-')
                cout << " ";
            else
                cout << endl;
        }
        cout << endl;
    }
#endif

    double starttime;


    IDWOperator idwOper(cellSize, idw_nbrPoints, idw_power, idw_buffer, blockSize);
    char* spatialrefWkt;

    idwOper.readSampleNums(inputFileName, &spatialrefWkt);
    vector<SamplePoint> samples;
    idwOper.readSamples(inputFileName, fldIdx, &spatialrefWkt, samples);
    idwOper.creatSampleBlocks(samples);
    vector<SamplePoint> vTemp(0);
    vTemp.swap(samples);

    RasterLayer<int> maskLayer("maskLayer");
    maskLayer.newLocalNbrhood();
    RasterLayer<double> idwLayer("idwLayer");
    idwLayer.newLocalNbrhood();

    if (decomposeBySapce) {//equal-area strategy
        idwOper.idwLayer(idwLayer, &spatialrefWkt);
        if(maskFileName){
            maskLayer.readFile(maskFileName, ROWWISE_DCMP);
            idwOper.maskLayer(maskLayer);
        }
        ComputeLayer<double> comptLayer("computeLayer");

        //if record execution time for preliminary experiment.
        //if true, the output file will be replaced by the time-recorded layer.
        if (writeLoadPath) {
            idwOper._writePreExpLoad = true;
        }
        if (myRank == 0) cout << "start computing" << endl;
        starttime = MPI_Wtime();
        idwOper.Run(); //fill idwLayer
    }
    else {//the proposed load-balancing strategy

         //Fill the spatial computational domain. This is a serial procedure.
        idwOper.idwLayerSerial(idwLayer, &spatialrefWkt);
        CoordBR subWorkBR;

        if (myRank == 0) cout << "start dcmp" << endl;
        if (readLoadPath) {//preliminary experiment mode, read the raster layer with recorded execution time.
            starttime = MPI_Wtime();
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.addRasterLayerSerial(&idwLayer);
            comptLayer.init(nullptr, granularity);
            comptLayer.readComputeLoadFile(readLoadPath);

            comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR); // Decompose the spatial computational domain.
        }
        else {//estimate function mode
            RasterLayer<int> fullMaskLayer("fullMaskLayer");
            if(maskFileName) {
                fullMaskLayer.newLocalNbrhood();
                fullMaskLayer.readGlobalFileSerial(maskFileName);
                idwOper.maskLayer(fullMaskLayer);
            }
            starttime = MPI_Wtime();
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.addRasterLayerSerial(&idwLayer);
            comptLayer.init(nullptr, granularity);
            IdwTransformation idwTrans(&comptLayer, &idwOper);
            idwTrans.run();
            comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR); // Decompose the spatial computational domain.
            if (writeLoadPath) {
                comptLayer.writeComputeIntensityFile(writeLoadPath);
            }
        }

        if (myRank == 0) cout << "dcmp time is " << MPI_Wtime() - starttime << endl;
        cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows() << endl;

        // Decompose the data domain using the extent of spatial computational subdomains
        idwOper.idwLayer(idwLayer, &spatialrefWkt, subWorkBR);

        if(maskFileName) {
            maskLayer.readFile(maskFileName, subWorkBR, ROWWISE_DCMP);
            idwOper.maskLayer(maskLayer);
        }

        if (myRank == 0) cout << "start computing" << endl;
        starttime = MPI_Wtime();
        idwOper.Run();
    }
    if (myRank == 0)
        cout << "compute time is " << MPI_Wtime() - starttime << endl << endl;
    idwLayer.writeFile(outputFileName);

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "write done." << endl;

    Application::END();
    return 0;
}
