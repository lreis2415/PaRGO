#ifndef SLOPEOPERATOR_H
#define SLOPEOPERATOR_H

#include "utility.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
#define Eps 0.0000001

class SlopeOperator : public RasterOperator<double> 
{
  public:
    SlopeOperator()
      :RasterOperator<double>(),
       cellSize(0), noData(-9999.), num(0),
       _pDEMLayer(nullptr), _pSlopeLayer(nullptr), _pDEMNbrhood(nullptr){}
   
    ~SlopeOperator() DEFAULT;
  
    void demLayer(RasterLayer<double> &layerD);
	void slopeLayer(RasterLayer<double> &layerD);

	bool isTermination() OVERRIDE;
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;

  protected:
	double cellSize;
	double noData;
	int num;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pSlopeLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif