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
#include "computLayer.h"
#include "application.h"
#include "idwOperator.h"
#include "communication.h"
#include "deComposition.h"

using namespace std;
using namespace GPRO;


int main(int argc, char *argv[]) 
{
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
	Application::START(MPI_Type, argc, argv); //init

	char *samplefilename, *outputfilename;
	char *dataNeighbor, *compuNeighbor;
	float cellSize;
	int fldIdx, idw_nbrPoints, idw_power, idw_buffer;	//暂时都定义为int可改为浮点型
	//int threadNUM;
	if( argc != 10 )
	{
		cout<<"Please Check the input paarameters!"<<endl;
		return 0;
	}
	samplefilename = argv[1];
	outputfilename = argv[2];
	dataNeighbor = argv[3];	//1*1邻域
	compuNeighbor = argv[4]; 
	cellSize = atof(argv[5]);	//待插值栅格分辨率
	fldIdx = atoi(argv[6]);	//矢量数据属性值所在列
	idw_power = atoi(argv[7]);	//反距离加权幂，通常取2
	idw_nbrPoints = atoi(argv[8]);	//搜索邻近点数
	idw_buffer = atof(argv[9]);	//最大搜索半径
	//threadNUM = atoi(argv[9]);

	//omp_set_num_threads(threadNUM);

	int myRank, process_nums;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	double starttime;
	double endtime;

	int blockGrain = 100;	//granularity,样点以块存放的粗网格粒度，以栅格分辨率为基本单位；用户根据数据指定
	IDWOperator idwOper(cellSize, idw_nbrPoints, idw_power, idw_buffer, blockGrain);
	char* spatialrefWkt;	//投影信息
	int sample_nums;
	sample_nums = idwOper.readSampleNums( samplefilename, &spatialrefWkt );	//获取样点数目，idwOper.sample_nums
	double **pAllSamples=(double **)malloc(sample_nums*sizeof(double *));//为采样点数组申请存储空间，存入子块后就释放
	for (int k=0; k<sample_nums; k++)
		pAllSamples[k]=(double *)malloc(3*sizeof(double));
	if (pAllSamples==NULL)
	{
		cout<<"Faliure memory request!"<<endl;
		return 0;
	}
	idwOper.readSamples( samplefilename, fldIdx, &spatialrefWkt, pAllSamples );	//读取样点，并更新了idwOper.glb_extent
	//可获取idwLayer的坐标范围,放入idwOper.sample_extent
	idwOper.creatSampleBlocks(pAllSamples);	//遍历pAllSamples，分块存入idwOper._pSampleBlocks成员
	//cout<<"creatSampleBlocks() done."<<endl;

	//以粗网格形式组织样点，数据成员行列数，每个栅格上是一系列样点
	RasterLayer<double> idwLayer("idwLayer");
	idwLayer.readNeighborhood(dataNeighbor);
	//equal row dcmp based on region
	idwOper.idwLayer(idwLayer, &spatialrefWkt);	//先将idwOperator的数据成员指向idwLayer图层，再借此创建idwLayer的基本元数据
	//cout<<"idwLayer metadata initialized."<<endl;
	//创建邻域类的临时对象，根据本图层的元数据直接划分,是否可行待定？
	starttime = MPI_Wtime();
	idwOper.Run();	//运行，结果写在idwLayer的cellspace中
	cout<<"idw compute done."<<endl;
	idwLayer.rowWriteFile(outputfilename, true);

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	if (myRank==0)
		cout<<"compute time is "<<endtime-starttime<<endl;

	//MPI_Barrier(MPI_COMM_WORLD);

	////balanced row dcmp based on compute burden
	//vInputLayers[0]->readNeighborhood(dataNeighbor);
	//int* pDcmpIdx = new int[process_nums*4];	//MPI不允许广播自定义类型，只能先写成角行列号
	//if( myRank==0 )	//这一步不应该交给用户，考虑内置到computLayer类去
	//{
	//	starttime = MPI_Wtime();
	//	vInputLayers[0]->readGlobalFile(vInputnames[0]);	//参与统计计算域的数据图层具有了完整的元数据,重复调用readFile时要保证clear并重新new
	//	vector<RasterLayer<double>* > inputLayers;
	//	inputLayers.push_back( vInputLayers[0] );		//若有多个图层参与计算域构建，这里push_back多次
	//	ComputLayer<double> comptLayer( inputLayers, "computLayer" );
	//	comptLayer.readNeighborhood(compuNeighbor);
	//	const int compuSize = 10;	//计算域图层分辨率是数据图层的10倍,粒度用户指定，这里暂定为10
	//	//获取负载均衡的划分，结果返回给vDcmpIdx,划分方式由第二个参数指定
	//	comptLayer.getCompuLoad( pDcmpIdx, ROWWISE_DCMP,compuSize, process_nums );	
	//	comptLayer.writeComptFile(outputfilename);
	//	endtime = MPI_Wtime();
	//	cout<<myRank<<" dcmp time is "<<endtime-starttime<<endl;
	//}
	////MPI_Barrier(MPI_COMM_WORLD);
	////依vDcmpIdx，主进程给各个进程广播其工作空间范围
	//MPI_Bcast(pDcmpIdx,process_nums*4,MPI_INT,0,MPI_COMM_WORLD);
	//CellCoord nwCorner(pDcmpIdx[4*myRank], pDcmpIdx[4*myRank+1]);
	//CellCoord seCorner(pDcmpIdx[4*myRank+2], pDcmpIdx[4*myRank+3]);
	//CoordBR subWorkBR(nwCorner, seCorner);
	//delete []pDcmpIdx;
	////MPI_Barrier(MPI_COMM_WORLD);
	//cout<<myRank<<" "<<subWorkBR<<" "<<subWorkBR.maxIRow() - subWorkBR.minIRow()<<endl;	//rowComptDcmp based on compt-burden has been ok.

	//vInputLayers[0]->readFile(vInputnames[0], subWorkBR, ROWWISE_DCMP);
	//for(int i=1; i<vInputnames.size(); i++){
	//	vInputLayers[i]->readNeighborhood(dataNeighbor);
	//	vInputLayers[i]->readFile(vInputnames[i], subWorkBR, ROWWISE_DCMP);
	//}
	//fcmLayer.copyLayerInfo(*vInputLayers[0]); //创建输出图层
	//for(int i=0; i<clusterNum; i++){
	//	vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
	//}

	cout<<"write done."<<endl;

	Application::END();
	//system("pause");
	return 0;
}
