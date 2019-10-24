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
#include "fcmOperator.h"
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

	//...
	char* inputfilenames;
	char* dataNeighbor;
	char* compuNeighbor;
	char* outputfilename;
	//int threadNUM;
	int clusterNum; //分类数目
	double maxIteration; //最大迭代次数
	double tolerance;//迭代阈值
	int wm;//加权指数
	if (argc>0 && argc < 10)
	{
		inputfilenames = argv[1];
		dataNeighbor = argv[2];
		compuNeighbor = argv[3]; 
		outputfilename = argv[4];
		clusterNum =atoi(argv[5]);
		maxIteration =atoi(argv[6]);
		tolerance = atof(argv[7]);
		wm = atof(argv[8]);		
		//threadNUM = atoi(argv[5]);
	}
	//omp_set_num_threads(threadNUM);

	int myRank, process_nums;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	double starttime;
	double endtime;

	//字符串解析输入文件名
	vector<char *> vInputnames;	//输入文件名，待解析
	vector<RasterLayer<double> *> vInputLayers;
	//int imageNum;	//输入的影像数目
	char* token = strtok(inputfilenames,",");
	char *filename;
	while(NULL != token)
	{
		filename=token;
		vInputnames.push_back(filename);	//名字留着后面读取用
		RasterLayer<double> *pLayer = new RasterLayer<double>("none");
		vInputLayers.push_back(pLayer);
		token=strtok(NULL,",");
	}
	if( vInputnames.empty() || clusterNum==0 || maxIteration==0 ){
		return 1;
	}
	RasterLayer<double> fcmLayer("fcmLayer");//创建分类输出图层fcmLayer
	//预定义分类图层
	char** pDegLayerName = new char*[clusterNum];
	//string degLayerName;
	vector<RasterLayer<double> *> vDegreeLayer;
	for(int i=0;i<clusterNum;i++)
	{
		pDegLayerName[i] = new char[50];
		sprintf(pDegLayerName[i],"degreeLayer%d.tif",i);	//输出用名字
		RasterLayer<double> *pLayer = new RasterLayer<double>(pDegLayerName[i]);
		vDegreeLayer.push_back(pLayer);
	}

	//equal row dcmp based on region
	for(int i=0; i<vInputnames.size(); i++){
		vInputLayers[i]->readNeighborhood(dataNeighbor);
		vInputLayers[i]->readFile(vInputnames[i], ROWWISE_DCMP);
	}
	fcmLayer.copyLayerInfo(*vInputLayers[0]); //创建输出图层
	for(int i=0; i<clusterNum; i++){
		vDegreeLayer[i]->copyLayerInfo(*vInputLayers[0]);
	}

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

	starttime = MPI_Wtime();
	FCMOperator fcmOper;
	fcmOper.initialization(vInputLayers.size(), clusterNum, maxIteration, tolerance, wm);
	fcmOper.inputLayer(vInputLayers);
	fcmOper.fcmLayer(fcmLayer);
	fcmOper.degLayer(vDegreeLayer);
	fcmOper.Run();

/*
	RasterLayer<double> demLayer("demLayer");
	RasterLayer<double> pitLayer("pitLayer"); //创建输出图层
	demLayer.readNeighborhood(dataNeighbor);
	
	//vector<CoordBR> vDcmpIdx;	//存放划分位置索引，主进程求得后广播给各进程;MPI不允许广播自定义类型，暂舍;
	//vDcmpIdx用什么数据类型合适，可在computLayer内转换格式，对用户透明，应尽量保证用户易用
	int* pDcmpIdx = new int[process_nums*4];	//MPI不允许广播自定义类型，只能先写成角行列号
	if( myRank==0 )	//这一步不应该交给用户，考虑内置到computLayer类去
	{
		starttime = MPI_Wtime();
		demLayer.readGlobalFile(inputfilename0);	//参与统计计算域的数据图层具有了完整的元数据,重复调用readFile时要保证clear并重新new
		vector<RasterLayer<double>* > inputLayers;
		inputLayers.push_back( &demLayer );		//若有多个图层参与计算域构建，这里push_back多次
		
		ComputLayer<double> comptLayer( inputLayers, "computLayer" );
		comptLayer.readNeighborhood(compuNeighbor);
		const int compuSize = 10;	//计算域图层分辨率是数据图层的10倍,粒度用户指定，这里暂定为10
		//获取负载均衡的划分，结果返回给vDcmpIdx,划分方式由第二个参数指定
		comptLayer.getCompuLoad( pDcmpIdx, ROWWISE_DCMP,compuSize, process_nums );	
		//for( int i=0; i<process_nums*4; i+=4 ){
		//	//it's ok here
		//	cout<<pDcmpIdx[i]<<" "<<pDcmpIdx[i+1]<<" "<<pDcmpIdx[i+2]<<" "<<pDcmpIdx[i+3]<<endl;
		//}
		//comptLayer.writeComptFile(outputfilename);	//可选,目前仅支持串行写

		endtime = MPI_Wtime();
		cout<<myRank<<" dcmp time is "<<endtime-starttime<<endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//依vDcmpIdx，主进程给各个进程广播其工作空间范围
	MPI_Bcast(pDcmpIdx,process_nums*4,MPI_INT,0,MPI_COMM_WORLD);
	CellCoord nwCorner(pDcmpIdx[4*myRank], pDcmpIdx[4*myRank+1]);
	CellCoord seCorner(pDcmpIdx[4*myRank+2], pDcmpIdx[4*myRank+3]);
	CoordBR subWorkBR(nwCorner, seCorner);
	delete []pDcmpIdx;
	MPI_Barrier(MPI_COMM_WORLD);
	cout<<myRank<<" "<<subWorkBR<<endl;	//rowComptDcmp based on compt-burden has been ok.

	demLayer.readFile(inputfilename0, subWorkBR, ROWWISE_DCMP);
	pitLayer.copyLayerInfo(demLayer);
*/
/*
	//method for equal decomposition based on region
	RasterLayer<double> demLayer("demLayer");
	RasterLayer<double> pitLayer("pitLayer"); //创建输出图层
	double starttime;
	double endtime;
	demLayer.readNeighborhood(dataNeighbor);
	demLayer.readFile(inputfilename0);
	pitLayer.copyLayerInfo(demLayer);	
*/

	//starttime = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	if (myRank==0)
		cout<<"compute time is "<<endtime-starttime<<endl;
	fcmLayer.writeFile(outputfilename);
	//for( size_t i = 0; i < vDegreeLayer.size(); ++i ){
	//	vDegreeLayer[i]->writeFile(pDegLayerName[i]);
	//}
	cout<<"write done."<<endl;

	Application::END();
	//system("pause");
	return 0;
}
