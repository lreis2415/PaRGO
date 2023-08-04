#ifndef LPITREMOVEOPERATOR_H
#define LPITREMOVEOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.000001

class PitRemoveOperator : public RasterOperator<double> 
{
  public:
    PitRemoveOperator()
      :RasterOperator<double>(),
       _pDEMLayer(0), _pwDEMLayer(0), _pwsLayer(0),_pwkLayer(0),_pOutLayer(0),_xSize(0), _ySize(0), num(0), flag(true){}
   
    ~PitRemoveOperator() {}

  
    void demLayer(RasterLayer<double> &layerD);
	void wdemLayer(RasterLayer<double> &layerD);
	void getarea(int &minrow,int &mincol,int &nrow,int &ncol,int _g,int buf,int id);
	void wslayer(RasterLayer<double> &layerD);
	void init(RasterLayer<double> &layerD);
	void writesca(RasterLayer<double> &layerD,CellCoord nw);
	void outLayer(RasterLayer<double>& layerD);
	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord,bool operFlag);

  public:
	double cellSize;
	double noData;
	int _xSize;	//rows in this processor
	int _ySize;
	int num;
	bool flag;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pwDEMLayer;
	Neighborhood<double> *_pDEMNbrhood;
	RasterLayer<double> *_pOutLayer;
	RasterLayer<double> *_pwsLayer;
	RasterLayer<double> *_pwkLayer;
};

#endif