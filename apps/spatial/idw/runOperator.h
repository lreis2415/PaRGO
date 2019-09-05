#ifndef RUNOPERATOR_H
#define RUNOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001
class RUNOperator : public RasterOperator<double> 
{
public:
	RUNOperator()
		:RasterOperator<double>(),
		_pDEMLayer(0), _pPitLayer(0), num(0), flag(true) {}
	//变量的初始化

	~RUNOperator() {}

	void demLayer(RasterLayer<double> &layerD);
	void pitLayer(RasterLayer<double> &layerD);
	
	virtual bool isTermination();

	virtual bool Operator(const CellCoord &coord, bool operFlag);


private:
	int _cellSize;
	int _xSize, _ySize;	//当前进程中DEM块的行数
	int _nRows, _nCols; //输入图层总行数和总列数
	double _noData;
	int _rank;
	int num;//迭代次数
	bool flag;

	RasterLayer<double> * _pDEMLayer;
	RasterLayer<double> * _pPitLayer;

	Neighborhood<double> *_pDEMNbrhood;

};

#endif