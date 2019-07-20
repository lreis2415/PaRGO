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

/*!
 * \enum SlopeAlgor
 */
enum SlopeAlgor {
    FD,  ///< Third-order finite difference weighted by reciprocal of squared distance
    FFD, ///< Frame finite difference
	MD,  ///< Maximum downslope
	SD,  ///< Simple difference
	SFD, ///< Second-order finite difference
	TFD, ///< Third-order finite difference
	TFDW ///< Third-order finite difference weighted by reciprocal of distance
};

SlopeAlgor GetSlopeAlgorithm(const string& arg);

class SlopeOperator : public RasterOperator<double> 
{
  public:
    SlopeOperator()
      :RasterOperator<double>(),
       cellSize(0), noData(-9999.), num(0), calcAlgor(FD),
       _pDEMLayer(nullptr), _pSlopeLayer(nullptr), _pDEMNbrhood(nullptr){}
   
    ~SlopeOperator() DEFAULT;
  
    void demLayer(RasterLayer<double> &layerD);
	void slopeLayer(RasterLayer<double> &layerD);
	void calcAlgorithm(SlopeAlgor algor);

	bool isTermination() OVERRIDE;
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;

  protected:
	double cellSize;
	double noData;
	int num;
	SlopeAlgor calcAlgor;
	RasterLayer<double> *_pDEMLayer;
	RasterLayer<double> *_pSlopeLayer;
	Neighborhood<double> *_pDEMNbrhood;
};

#endif