#ifndef TWIOPERATOR_H
#define TWIOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class TWIOperator : public RasterOperator<double> 
{
  public:
    TWIOperator()
      :RasterOperator<double>(),
	_pSCALayer(0), _pSlopeLayer(0),_pTWILayer(0), flag(true){}
  
    ~TWIOperator() {}

	void SCALayer(RasterLayer<double> &layerD);
	void slopeLayer(RasterLayer<double> &layerD);
	void twiLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();

    virtual bool Operator(const CellCoord &coord, bool Operflag);

  public:
	double _cellSize;
	double _noData;
	bool flag;
	int _iNeighborCells;

	RasterLayer<double> *_pSCALayer;
	RasterLayer<double> *_pSlopeLayer;
	RasterLayer<double> *_pTWILayer;
	Neighborhood<double> *_pNbrhood;
};

#endif