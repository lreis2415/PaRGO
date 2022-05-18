#ifndef DINFOPERATOR_H
#define DINFOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class DinfOperator : public RasterOperator<double> 
{
  public:
    DinfOperator()
      :RasterOperator<double>(),
	  cellSize(0), noData(-9999.), num(0),_maxRow(0),_maxCol(0),
       _pDinfLayer(nullptr), _pSCALayer(nullptr), _pDinfNbrhood(nullptr),_degreeLayer("degreeLayer"){}
 
    ~DinfOperator() {}

  
    void dinfLayer(RasterLayer<double> &layerD);
	void scaLayer(RasterLayer<double> &layerD);
	


	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);
	//virtual bool Operator(const CellCoord& coord, bool operFlag);
  protected:
	double cellSize;
	double noData;
	int num;
	int _maxRow;
	int _maxCol;

	RasterLayer<double> _degreeLayer;
	RasterLayer<double> *_pDinfLayer;
	RasterLayer<double> *_pSCALayer;
	Neighborhood<double> *_pDinfNbrhood;
};

#endif