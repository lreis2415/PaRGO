#ifndef STREAMOPERATOR_H
#define STREAMOPERATOR_H

#include "utility.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;
#define Eps 0.0000001

class streamOperator : public RasterOperator<double> 
{
  public:
    streamOperator()
      :RasterOperator<double>(),
       _noData(-9999.), num(0),maxval(0.),length(0.),runtime(0.),_maxrow(0),_maxcol(0),_contribs("contribs"),_porderout(0),
       _pDirLayer(0), _pnetLayer(0),s1(0.),s2(0.),s1sq(0.),s2sq(0.),n1(0),n2(0),totalarea(0.){}
   
	~streamOperator() {};
  
    //void demLayer(RasterLayer<double> &layerD);
	void netLayer(RasterLayer<double> &layerD);
	void dirLayer(RasterLayer<double> &layerD);
	//void areaLayer(RasterLayer<double> &layerD);
	void orderLayer(RasterLayer<double> &layerD);
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
	//double &opt;
	double runtime;

	RasterLayer<double> *_pnetLayer;
	RasterLayer<double> *_pDirLayer;
	RasterLayer<double> _contribs;
	RasterLayer<double> *_porderout;
	Neighborhood<double> *_pnetNbrhood;
	//Neighborhood<double> *_pDirNbrhood;
};

#endif