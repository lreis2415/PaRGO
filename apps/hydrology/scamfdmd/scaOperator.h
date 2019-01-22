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
       _pSCALayer(0), _pNbrhood(0), num(0), _maxRow(0), _maxCol(0), _iNeighborCells(0),
	   _degreeLayer("degreeLayer") {}
 
    ~SCAOperator() {}

  
    void mfdLayer(vector<RasterLayer<double> *> &layerV);
	void scaLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  protected:
	double _cellSize;
	double _noData;
	int num;
	int _maxRow;
	int _maxCol;
	int _iNeighborCells;

	RasterLayer<double> _degreeLayer;
	vector<RasterLayer<double> *> _weightLayerVec;
	vector<CellSpace<double>* > weightLs;
	RasterLayer<double> *_pSCALayer;
	Neighborhood<double> *_pNbrhood;
};

#endif