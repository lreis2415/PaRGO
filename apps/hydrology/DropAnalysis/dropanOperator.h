#ifndef DROPANOPERATOR_H
#define DROPANOPERATOR_H

#include "utility.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
#define Eps 0.0000001

class DropanOperator : public RasterOperator<double> 
{
  public:
    DropanOperator(double th,double &opt)
      :RasterOperator<double>(),
       _noData(-9999.), num(0),maxval(0.),length(0.),_maxrow(0),_maxcol(0),_elevOut("elevOut"),_contribs("contribs"),_orderout("orderout"),
       _pDirLayer(0), _pssaLayer(0),_pDEMLayer(0),thresh(th),s1(0.),s2(0.),s1sq(0.),s2sq(0.),n1(0),n2(0),totalarea(0.),opt(opt),
		_pareaLayer(0){}
   
	~DropanOperator() {};
  
    void demLayer(RasterLayer<double> &layerD);
	void ssaLayer(RasterLayer<double> &layerD);
	void dirLayer(RasterLayer<double> &layerD);
	void areaLayer(RasterLayer<double> &layerD);

	virtual bool isTermination();
    virtual bool Operator(const CellCoord &coord, bool operFlag);
	
	//double maxval,totalarea;

  protected:
	double _noData;
	double _cellSize;
	int _rank;
	int num;
	int _maxrow;
	int _maxcol;
	double maxval;
	double length;
	double thresh;
	double s1,s2,s1sq,s2sq;
	long n1,n2;
	double totalarea;
	double &opt;

	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pssaLayer;
	RasterLayer<double> *_pDirLayer;
	RasterLayer<double> *_pareaLayer;
	RasterLayer<double> _elevOut;
	RasterLayer<double> _contribs;
	RasterLayer<double> _orderout;
	Neighborhood<double> *_pssaNbrhood;
	//Neighborhood<double> *_pDirNbrhood;
};

#endif