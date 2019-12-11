#ifndef RUNOPERATOR_H
#define RUNOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include "transformation.h"
#include "utility.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001

class FCMOperator : public RasterOperator<double> 
{
public:
	FCMOperator()
		:RasterOperator<double>(),
		 _iterNum(0), flag(true),
		block(0),nodataNums(0), dFlag(true),subval(0.0),totval(0.0),oldtval(0.0),partitionCoef(0.0),entropy(0),totpartitionCoef(0),totentropy(0)
		{}
	//变量的初始化

	~FCMOperator();

	void initialization(int iNum, int cNum, int maxIter,double toler, double m);
	void inputLayer(vector<RasterLayer<double> *> layerD);
	void fcmLayer(RasterLayer<double> &layerD);
	void degLayer(vector<RasterLayer<double> *> layerD);
	
	virtual bool isTermination();

	void createRandomIdx( int nums, int range, int* randomIdx );	//在range范围内产生nums个随机数，randomIdx返回
	void fnDistance(int curRow, int curCol, double* pInputVal); //计算距离
	void InitDegree(int curRow, int curCol);//计算隶属度
	void initRandomClusterCenters(double *clusterCenters);
	void assignMaxMembershipDegrees();
	virtual bool Operator(const CellCoord &coord, bool operFlag);


private:
	int _cellSize;
	int _xSize, _ySize;	//当前进程中DEM块的行数
	int _nRows, _nCols; //输入图层总行数和总列数
	double _noData;
	int _rank;
	int _iterNum;//迭代次数
	bool flag;
	double starttime, endtime;
	double tmpSumTime1, tmpSumTime2;
protected:
	int clusterNum; //分类数目
	int maxIteration; //最大迭代次数
	double tolerance;//迭代阈值
	double weight;//加权指数
	int imageNum;//输入的影像数目

	int block;	//各进程分块数据量
	int nodataNums;
	bool dFlag;		//用来标记是否第一次迭代
	double* centerVal;	//聚类中心属性值，一维数组
	int* centerIndex;	//聚类中心位置索引 clusterNum * 1大小数组，xy坐标映射的一维
	double*** degree;	//隶属度数组，三维
	double*** dist;	//存放各个点与各类中心的距离
	double subval;	//存放本次迭代的约束函数值
	double totval;	//主进程汇总各进程，得到阈值和
	double oldtval;	//上次迭代的阈值和
	double *sumNumerator;	//分子求和
	double *sumDenominator;		//分母求和
	double *totNumerator;	//分子归并
	double *totDenominator;		//分母归并
	double partitionCoef;	//分割系数
	double entropy;	//熵
	double totpartitionCoef;	//规约后的分割系数
	double totentropy;	//规约后的熵

	vector<RasterLayer<double> *> _vInputLayer;
	RasterLayer<double> *_pFCMLayer;
	vector<RasterLayer<double> *> _vDegLayer;

	Neighborhood<double> *_pDEMNbrhood;

};

#endif