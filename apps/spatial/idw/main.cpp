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
        << "-dataNbr <data layer neighbor file> " << endl
        << "-computeNbr <compute layer neighbor file>" << endl
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


    //cout << "Example.1. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif" << endl;
    //cout << "Example.2. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif -mtd SD" << endl;
    //cout << "Example.3. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif" << endl;
    //cout << "Example.4. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif TFD" << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName;
    char* outputFileName;
    char* dataNeighbor;
    char* compuNeighbor;
    float cellSize;
    int fldIdx;
    int idw_nbrPoints;
    int idw_power;
    double idw_buffer; //暂时都定义为int可改为浮点型
    double blockSize; //或称granularity,样点以块存放的粗网格粒度，以米为单位
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
    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);
    if (myRank == 0) {
        cout<<"PaRGO-IDW. "<<process_nums<<" core(s)"<<endl;
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


    IDWOperator idwOper(cellSize, idw_nbrPoints, idw_power, idw_buffer, blockSize);
    char* spatialrefWkt; //投影信息
    int sample_nums;
    sample_nums = idwOper.readSampleNums(inputFileName, &spatialrefWkt); //获取样点数目，idwOper.sample_nums
    double** pAllSamples = (double **)malloc(sample_nums * sizeof(double *)); //为采样点数组申请存储空间，存入子块后就释放
    for (int k = 0; k < sample_nums; k++)
        pAllSamples[k] = (double *)malloc(3 * sizeof(double));
    if (pAllSamples == NULL) {
        cout << "Faliure memory request!" << endl;
        return 0;
    }
    idwOper.readSamples(inputFileName, fldIdx, &spatialrefWkt, pAllSamples); //读取样点，并更新了idwOper.glb_extent
    //可获取idwLayer的坐标范围,放入idwOper.sample_extent
    idwOper.creatSampleBlocks(pAllSamples); //遍历pAllSamples，分块存入idwOper._pSampleBlocks成员
    delete []pAllSamples;
    //以粗网格形式组织样点，数据成员行列数，每个栅格上是一系列样点
    RasterLayer<double> idwLayer("idwLayer");
    idwLayer.readNeighborhood(dataNeighbor);
    
    if (decomposeBySapce) {
        idwOper.idwLayer(idwLayer, &spatialrefWkt);
        ComputeLayer<double> comptLayer("computeLayer");
        if (writeLoadPath) {
            comptLayer.initSerial(&idwLayer, compuNeighbor, 1);
            idwOper.comptLayer(comptLayer);
        }
        starttime = MPI_Wtime();
        idwOper.Run(); //fill idwLayer
        if (writeLoadPath)
            comptLayer.writeComputeIntensityFileSerial(writeLoadPath);
    }
    else {
        idwOper.idwLayerSerial(idwLayer, &spatialrefWkt);
        const int compuSize = 10; //ComputeLayerGrain=10*DataLayerGrain (10 temporarily)
        CoordBR subWorkBR;
        if (readLoadPath) {
            ComputeLayer<double> comptLayer("computLayer");
            comptLayer.initSerial(&idwLayer, compuNeighbor, compuSize);
            comptLayer.readComputeLoadFile(readLoadPath);
            comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
        }
        else {
            ComputeLayer<double> comptLayer("computLayer");
            starttime = MPI_Wtime();

            comptLayer.initSerial(&idwLayer, compuNeighbor, compuSize);
            IdwTransformation idwTrans(&comptLayer, &idwOper);
            idwTrans.run();
            comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);

            endtime = MPI_Wtime();
            cout << myRank << " dcmp time is " << endtime - starttime << endl;

            if (writeLoadPath) {
                comptLayer.writeComputeIntensityFileSerial(writeLoadPath);
            }
        }

        cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows() << endl;
        idwOper.idwLayer(idwLayer, &spatialrefWkt, subWorkBR);
        starttime = MPI_Wtime();
        
        idwOper.Run();
    }
    endtime = MPI_Wtime();
    if (myRank == 0)
        cout << "compute time is " << endtime - starttime << endl << endl;
    idwLayer.writeFile(outputFileName);

    MPI_Barrier(MPI_COMM_WORLD);
    //cout << "write done." << endl;

    Application::END();
    return 0;
}
