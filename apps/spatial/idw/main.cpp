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

    cout << " Or use the Simple Usage (full-text usage is recommended): idw <sample points file> "
        "<data layer neighbor file> "
        "<compute layer neighbor file>"
        "<output raster file>"
        "<output resolution (length of one cell)(meter)>"
        "<index of the field to be calculated>"
        "<exponent of inverse distance>"
        "<maximum of neighbor points to search>"
        "<maximum of range to search>"
        "<domain decomposition method>"
        "[<path to write load file>]"
        "[<path to read load file>]"<< endl << endl;

    //cout << "Example.1. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif" << endl;
    //cout << "Example.2. slope -elev /path/to/elev.tif -nbr /path/to/moore.nbr -slp /path/to/slp.tif -mtd SD" << endl;
    //cout << "Example.3. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif" << endl;
    //cout << "Example.4. slope /path/to/elev.tif /path/to/moore.nbr /path/to/slp.tif TFD" << endl;

    exit(1);
}

int main(int argc, char *argv[]) 
{
	char *inputFileName;
	char *outputFileName;
	char *dataNeighbor;
	char *compuNeighbor;
	float cellSize;
	int fldIdx;
	int idw_nbrPoints;
	int idw_power;
	int idw_buffer;	//暂时都定义为int可改为浮点型
    bool decomposeBySapce; /// decomp by compute load if false
    char* writeLoadPath=nullptr;
    char* readLoadPath=nullptr;

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
                Usage("No argument followed '-computeNbr'!");
            }
        }
		else if (strcmp(argv[i], "-resolution") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                cellSize = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-computeNbr'!");
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
                Usage("No argument followed '-computeNbr'!");
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
                Usage("No argument followed '-computeNbr'!");
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
                Usage("No argument followed '-computeNbr'!");
            }
        }
		else if (strcmp(argv[i], "-searchRange") == 0) {
            simpleusage = false;
            i++;
            if (argc > i) {
                idw_buffer = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-computeNbr'!");
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
            dataNeighbor = argv[2];//1*1邻域
            compuNeighbor = argv[3];
            outputFileName = argv[4];
            cellSize = atoi(argv[5]);//待插值栅格分辨率
            fldIdx = atoi(argv[6]);//矢量数据属性值所在列
            idw_power = atoi(argv[7]);//反距离加权幂，通常取2
            idw_nbrPoints = atoi(argv[8]);//搜索邻近点数
            idw_buffer = atoi(argv[8]);//最大搜索半径
			if(argc>=10){
            	decomposeBySapce = ifDecomposeBySpace(argv[9]);
			}
            if(argc>=11){
                writeLoadPath = argv[10];
                readLoadPath = argv[10];
            }
            break;
        }
	}

	Application::START(MPI_Type, argc, argv);

	int myRank, process_nums;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name,&name_len);
    if(myRank==0)
        cout<<"pargo-idw. decompose by space: "<<decomposeBySapce<<". "<<process_nums<<" core(s)."<<endl;
    cout << "process " << myRank << " on " << processor_name << endl;
    
	double starttime;
	double endtime;

	int blockGrain = 5;	//granularity,样点以块存放的粗网格粒度，以栅格分辨率为基本单位；用户根据数据指定
	IDWOperator idwOper(cellSize, idw_nbrPoints, idw_power, idw_buffer, blockGrain);
	char* spatialrefWkt;	//投影信息
	int sample_nums;
	sample_nums = idwOper.readSampleNums( inputFileName, &spatialrefWkt );	//获取样点数目，idwOper.sample_nums
	double **pAllSamples=(double **)malloc(sample_nums*sizeof(double *));//为采样点数组申请存储空间，存入子块后就释放
	for (int k=0; k<sample_nums; k++)
		pAllSamples[k]=(double *)malloc(3*sizeof(double));
	if (pAllSamples==NULL)
	{
		cout<<"Faliure memory request!"<<endl;
		return 0;
	}
	idwOper.readSamples( inputFileName, fldIdx, &spatialrefWkt, pAllSamples );	//读取样点，并更新了idwOper.glb_extent
	//可获取idwLayer的坐标范围,放入idwOper.sample_extent
	idwOper.creatSampleBlocks(pAllSamples);	//遍历pAllSamples，分块存入idwOper._pSampleBlocks成员

	//以粗网格形式组织样点，数据成员行列数，每个栅格上是一系列样点
	RasterLayer<double> idwLayer("idwLayer");
	idwLayer.readNeighborhood(dataNeighbor);

	if(decomposeBySapce){
		//equal row dcmp based on region
		idwOper.idwLayer(idwLayer, &spatialrefWkt);	//先将idwOperator的数据成员指向idwLayer图层，再借此创建idwLayer的基本元数据
		//创建邻域类的临时对象，根据本图层的元数据直接划分,是否可行待定？
		ComputeLayer<double> comptLayer("computeLayer");
        if(writeLoadPath) {
            comptLayer.init(&idwLayer,compuNeighbor,10);
            //comptLayer.copyLayerInfo(idwLayer);
            //comptLayer.newMetaData(10);
            idwOper.comptLayer(comptLayer);
        }
		starttime = MPI_Wtime();
		idwOper.Run();	//运行，结果写在idwLayer的cellspace中
		if(writeLoadPath)
			comptLayer.writeFile(writeLoadPath); //测试用，写出捕捉到的计算时间
	}
	else{
		//balanced row dcmp based on compute burden
        const int compuSize = 10;	//计算域图层分辨率是数据图层的10倍,粒度用户指定，这里暂定为10
		idwOper.idwLayer(idwLayer, &spatialrefWkt);
		idwLayer._pMetaData->_domDcmpType = NON_DCMP;	//wyj 2019-11-12:暂时写在外面
		CoordBR subWorkBR;

        RasterLayer<double> idwGlobalLayer("globalLayer");
	    idwGlobalLayer.readNeighborhood(dataNeighbor);
     //   idwOper.idwLayer(idwGlobalLayer, &spatialrefWkt);
        idwGlobalLayer.copyLayerInfo(idwLayer);
        if(myRank==0) {
            idwGlobalLayer.newCellSpace(idwGlobalLayer.metaData()->_glbDims,0);            
        }
        idwGlobalLayer._pMetaData->_localworkBR.seCorner(idwGlobalLayer._pMetaData->_glbDims.nRows()-idwGlobalLayer.nbrhood()->nRows(),
            idwGlobalLayer._pMetaData->_glbDims.nCols()-idwGlobalLayer.nbrhood()->nCols());
        ComputeLayer<double> comptLayer(&idwGlobalLayer,compuSize,"computLayer");
        //comptLayer.setComputGrain(compuSize);
        comptLayer.copyLayerMetadata(idwGlobalLayer);
        comptLayer.readNeighborhood(compuNeighbor);
        if(readLoadPath) {
            if(myRank==0) {
                comptLayer.readFile(readLoadPath);  
            }else {
                MPI_Barrier( MPI_COMM_WORLD );
            }
            comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
        }else{
            if( myRank==0 )
            {
                starttime = MPI_Wtime();
                comptLayer.init(compuSize);
                IdwTransformation idwTrans(&comptLayer,&idwOper);
                idwTrans.run();
                comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);	
                if (writeLoadPath) {
                    comptLayer.writeComputeFile(writeLoadPath);
                }
                endtime = MPI_Wtime();
                cout<<myRank<<" dcmp time is "<<endtime-starttime<<endl;
            }else{
                ComputeLayer<double> comptLayer("untitled");
                comptLayer.getCompuLoad( ROWWISE_DCMP, process_nums, subWorkBR );
                if(writeLoadPath) {
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
        }
	    cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows()<< endl;
        idwLayer.layerDcmp(subWorkBR);
		idwOper.idwLayer(idwLayer, &spatialrefWkt,subWorkBR);
		starttime = MPI_Wtime();

		idwOper.Run();
	}
	endtime = MPI_Wtime();
	if (myRank==0)
		cout<<"compute time is "<<endtime-starttime<<endl<<endl;
	idwLayer.writeFile(outputFileName);
    
	MPI_Barrier(MPI_COMM_WORLD);
  //  if(myRank==2) {
  //      double t1,t2;
		//t1 = MPI_Wtime();
		//t2 = MPI_Wtime();
  //      while(t2-t1<600) {
  //          t2=MPI_Wtime();
  //      }
  //  }
	cout<<"write done."<<endl;

	Application::END();
	//system("pause");
	return 0;
}
