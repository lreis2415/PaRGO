#ifndef SCAOPERATOR_H
#define SCAOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>
#include <vector>
#include<omp.h>

#include "ogrsf_frmts.h"

#define Eps 1e-9
using namespace GPRO;

class SCAOperator : public RasterOperator<double>
{
public:
	SCAOperator()
		:RasterOperator<double>(),
		_pDEMLayer(0), _pSCALayer(0),/*_pTMPLayer(0), */num(0), flag(true){}

	~SCAOperator() {}


	void demLayer(RasterLayer<double> &layerD, int kc_meth2, float StepRatio2, int threadNUM2);
	void SCALayer(RasterLayer<double> &layerD);
	//void tmpLayer(RasterLayer<double> &layerD);
	virtual bool isTermination();

	inline double q1(double t){
		return 2.0*t*t / 3.0;
	}
	inline double q2(double t){
		return ((1 - t)*4.0*t + 2.0) / 3.0;
	}
	inline double q3(double t){
		return 2.0*(1 - t)*(1 - t) / 3.0;
	}

	virtual bool Operator(const CellCoord &coord);

public:
	double cellSize;
	double noData;
	int num;	//控制当前迭代次数
	bool flag;
	int kc_Meth;
	float StepR;
	int threadNum;
	double** zq;
	double** qb;
	double** diff;
	double* relDiffMin;
	double* relDiffMax;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pSCALayer;
	//RasterLayer<double> *_pTMPLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif