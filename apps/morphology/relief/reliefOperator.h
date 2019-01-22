#ifndef RELIEFOPERATOR_H
#define RELIEFOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
class reliefOperator : public RasterOperator<double> 
{
  public:
    reliefOperator()
      :RasterOperator<double>(),
       _pDEMLayer(0), _pReliefLayer(0), num(0){}
   
    ~reliefOperator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void reliefLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();

    virtual bool Operator(const CellCoord &coord);

  protected:
	int cellSize;
	int noData;
	int num;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pReliefLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif