#ifndef IDWOPERATOR_H
#define IDWOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001
#define  NODATA_DEFINE -9999

struct Sample_Point{
	int x;
	int y;
	double value;
};

struct extent_info {
	double minX;
	double maxX;
	double minY;
	double maxY;
};

struct Sample_block{
	//CoordBR _MBR;	//是否需要待定
	//extent_info blockExtent;
	vector<Sample_Point> sample_Points;
};

class IDWOperator : public RasterOperator<double> 
{
public:
	IDWOperator()
		:RasterOperator<double>(),
		 _iterNum(0), flag(true), _sample_nums(0)
		{}

	//变量的初始化
	IDWOperator( float cellsize, int nbrPoints, int idw_power, int bufferSize, int grain )
		:RasterOperator<double>(),
		_iterNum(0), flag(true), _sample_nums(0), _cellSize(cellsize), _nbrPoints(nbrPoints), _idw_power(idw_power), _idw_buffer(bufferSize), _blockGrain(grain), _noData(NODATA_DEFINE)
	{}

	~IDWOperator();

	int readSampleNums( const char* filename, char** pSpatialRefWkt );
	bool readSamples( const char* filename, int fieldIdx, char** pSpatialRefWkt, double **Sample_Array );
	void creatSampleBlocks( double **pSamples );

	void idwLayer(RasterLayer<double> &layerD, char** pSpatialRefWkt);
	virtual bool isTermination();

	int searchNbrSamples( const int subMinRow, int cellRow, int cellCol, double *nbrSamples);

	virtual bool Operator(const CellCoord &coord, bool operFlag);


private:
	float _cellSize;
	int _xSize, _ySize;	//当前进程中DEM块的行数
	int _nRows, _nCols; //输入图层总行数和总列数
	double _noData;
	int _myRank;
	int _iterNum;//迭代次数
	bool flag;
	double starttime, endtime;

protected:
	int _nbrPoints;	//插值邻域样点数
	int _idw_power;
	int _idw_buffer;
	extent_info _glb_extent;	//所构建的插值栅格图层全区范围，根据样点数据范围外扩确定
	extent_info _sub_extent;	//本进程栅格范围
	int _sample_nums;	//全区总样点数量
	int _blockGrain;
	Sample_block* _pSampleBlocks;	//以粗网格块,按行组织的样点;每个进程都存了全部块，析构函数中要释放；

	RasterLayer<double> *_pIDWLayer;

};

#endif