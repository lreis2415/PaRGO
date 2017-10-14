#ifndef MFDOPERATOR_H
#define MFDOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class MFDOperator : public RasterOperator<double> 
{
  public:
    MFDOperator()
      :RasterOperator<double>(),
       _pDEMLayer(0), _pDEMNbrhood(0), slpExp(0.0), num(0){}
	MFDOperator( double pExp )
		:RasterOperator<double>(),
		_pDEMLayer(0), _pDEMNbrhood(0), slpExp(pExp), num(0){}
  
    ~MFDOperator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void weightLayers(vector<RasterLayer<double> *> &layers);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  protected:
	double _cellSize;
	double _noData;
	double slpExp;
	int num;
	RasterLayer<double> *_pDEMLayer;
	vector<RasterLayer<double> *> _weightLayers;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif