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
	  cellSize(0), noData(-9999.), num(0),
       _pDEMLayer(nullptr), _pslpLayer(nullptr),_pDinfLayer(nullptr), _pDEMNbrhood(nullptr){}
 
    ~DinfOperator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void scaLayer(RasterLayer<double> &layerD);
	void slpLayer(RasterLayer<double> &layerD);
	void dinfLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  protected:
	double cellSize;
	double noData;
	int num;
	//int _maxRow;
	//int _maxCol;

	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pslpLayer;
	RasterLayer<double> *_pDinfLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif