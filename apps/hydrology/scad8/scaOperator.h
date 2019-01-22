#ifndef SCAOPERATOR_H
#define SCAOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class SCAOperator : public RasterOperator<double> 
{
  public:
    SCAOperator()
      :RasterOperator<double>(),
       _pD8Layer(0), _pSCALayer(0), num(0), _maxRow(0), _maxCol(0),
	   _degreeLayer("degreeLayer") {}
 
    ~SCAOperator() {}

  
    void d8Layer(RasterLayer<double> &layerD);
	void scaLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  protected:
	double _cellSize;
	double _noData;
	int num;
	int _maxRow;
	int _maxCol;

	RasterLayer<double> _degreeLayer;
	RasterLayer<double> *_pD8Layer;
	RasterLayer<double> *_pSCALayer;
	Neighborhood<double> *_pD8Nbrhood;
};

#endif