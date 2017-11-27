#ifndef D8OPERATOR_H
#define D8OPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class D8Operator : public RasterOperator<double> 
{
  public:
    D8Operator()
      :RasterOperator<double>(),
       _pDEMLayer(0), _pD8Layer(0), dirType(0), num(0){}
	D8Operator( int type )
		:RasterOperator<double>(),
		_pDEMLayer(0), _pD8Layer(0), dirType(type), num(0){}
  
    ~D8Operator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void d8Layer(RasterLayer<double> &layerD);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  protected:
	double cellSize;
	double noData;
	int dirType;
	int num;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pD8Layer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif