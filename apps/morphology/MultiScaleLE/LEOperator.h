#ifndef LEOPERATOR_H
#define LEOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class LEOperator : public RasterOperator<double> 
{
  public:
    LEOperator()
      :RasterOperator<double>(),
       _pDEMLayer(0), _pLELayer(0), num(0){}
   
    ~LEOperator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void LELayer(RasterLayer<double> &layerD);
	virtual bool isTermination();

    bool Operator(const CellCoord &coord,bool operFlag);

    int GetRowNum() { return sqrt((double)(_pDEMNbrhood->size() - 1) / 2); }

    int GetCurrentScale() { return _icurrentScale; };
    void SetCurrentScale(int scale) { _icurrentScale = scale; }
  protected:
	int cellSize;
	int noData;
	int num;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pLELayer;
	Neighborhood<double> *_pDEMNbrhood;
private:
    int _icurrentScale;

};

#endif