#ifndef RESAMPLEOPERATOR_H
#define RESAMPLEOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class ResampleOperator : public RasterOperator<double> 
{
public:
	ResampleOperator()
		:RasterOperator<double>(){};


    void setInputLayer(RasterLayer<double>& inputLayer);
    void setOutputLayer(RasterLayer<double>& outputLayer);
    void setG(int g){_g=g;}
private:
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;
    
	RasterLayer<double> *_inputLayer;
	RasterLayer<double> *_outputLayer;
    int _g;//granularity
    


};

#endif