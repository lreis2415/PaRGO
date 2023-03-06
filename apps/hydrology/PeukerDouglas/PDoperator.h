#ifndef PEUKERDOUGLAS_H
#define PEUKERDOUGLAS_H

#include "utility.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
#define Eps 0.0000001

class PDOperator : public RasterOperator<double> 
{
  public:
    PDOperator()
      :RasterOperator<double>(),
       noData(-9999.), num(0),_maxRow(0), _maxCol(0),_smoLayer("smoLayer"),
       _pDEMLayer(nullptr), _pucgLayer(nullptr),_pDEMNbrhood(nullptr){}
   
    ~PDOperator() DEFAULT;
  
    void demLayer(RasterLayer<double> &layerD);
	void ucgLayer(RasterLayer<double> &layerD);

	bool isTermination() OVERRIDE;
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;

  protected:
	double noData;
	int num;
	int _maxRow;
	int _maxCol;

	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pucgLayer;//upward curved grid
	RasterLayer<double> _smoLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif