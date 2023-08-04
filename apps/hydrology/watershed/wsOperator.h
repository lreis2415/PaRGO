#ifndef WSOPERATOR_H
#define WSOPERATOR_H

#include "utility.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
#define Eps 0.0000001

class wsOperator : public RasterOperator<double> 
{
  public:
    wsOperator()
      :RasterOperator<double>(),
       _noData(-9999.), num(0),maxval(0.),length(0.),runtime(0.),_maxrow(0),_maxcol(0),_contribs("contribs"),_pWatershed(0),myrank(0),ranksize(0),
       _pDirLayer(0), _psmLayer(0),n1(0),n2(0){}
   
	~wsOperator() {};
  
    //void demLayer(RasterLayer<double> &layerD);
	void smLayer(RasterLayer<double> &layerD);
	void dirLayer(RasterLayer<double> &layerD);
	//void areaLayer(RasterLayer<double> &layerD);
	void wsLayer(RasterLayer<double> &layerD);
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
	int myrank;
	int ranksize;
	double maxval;
	double length;
	double thresh;
	double s1,s2,s1sq,s2sq;
	long n1,n2;
	double totalarea;
	//double &opt;
	double runtime;

	RasterLayer<double> *_psmLayer;
	RasterLayer<double> *_pDirLayer;
	RasterLayer<double> _contribs;
	RasterLayer<double> *_pWatershed;
	Neighborhood<double> *_psmNbrhood;
	//Neighborhood<double> *_pDirNbrhood;
};

#endif