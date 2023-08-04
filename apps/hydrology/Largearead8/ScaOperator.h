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
       _pD8Layer(0), _pSCALayer(0), _pwsLayer(0),_pwkLayer(0),_pOutLayer(0),num(0), _maxRow(0), _maxCol(0),
	   _degreeLayer("degreeLayer") {}
 
    ~SCAOperator() {}

	void getarea(int &minrow,int &mincol,int &nrow,int &ncol,int _g,int buf,int id);
    void d8Layer(RasterLayer<double> &layerD);
	void scaLayer(RasterLayer<double> &layerD);
	void wslayer(RasterLayer<double> &layerD);
	void initsca(RasterLayer<double> &layerD);
	void writesca(RasterLayer<double> &layerD,CellCoord nw);
	void outLayer(RasterLayer<double>& layerD);
	//void initsca(); 

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
	RasterLayer<double> *_pOutLayer;
	RasterLayer<double> *_pwsLayer;
	RasterLayer<double> *_pwkLayer;
	RasterLayer<double> *_pSCALayer;
	Neighborhood<double> *_pD8Nbrhood;
};

#endif