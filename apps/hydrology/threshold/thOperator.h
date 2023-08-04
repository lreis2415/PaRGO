#ifndef THOPERATOR_H
#define THOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.000001

class thOperator : public RasterOperator<double> 
{
  public:
    thOperator(double th)
      :RasterOperator<double>(),
       _pSCALayer(0), _pnetLayer(0), _xSize(0), _ySize(0), num(0), flag(true),th(th){}
   
    ~thOperator() {}

  
    void scaLayer(RasterLayer<double> &layerD);
	void netLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  public:
	double cellSize;
	double noData;
	int _xSize;	//rows in this processor
	int _ySize;
	int num;
	bool flag;
	double th;
	RasterLayer<double> *_pSCALayer;
	RasterLayer<double> *_pnetLayer;
	Neighborhood<double> *_pSCANbrhood;
};

#endif