#ifndef DEMO1_H
#define DEMO1_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace GPRO;

class ReclassifyOperator : public RasterOperator<double> 
{
public:
	ReclassifyOperator()
		:RasterOperator<double>(){};

    void setInputLayer(RasterLayer<double>& inputLayer);
    void setOutputLayer(RasterLayer<double>& outputLayer){_outputLayer=&outputLayer;}
    void setLevels(vector<double>* dLevels){levels=dLevels;}
    vector<double>* levels;

private:
    bool Operator(const CellCoord &coord, bool operFlag) OVERRIDE;

	RasterLayer<double> *_inputLayer;
	RasterLayer<double> *_outputLayer;
    


};

#endif