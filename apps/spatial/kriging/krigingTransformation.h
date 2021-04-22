#ifndef KRIGINGTRANSFORMATION_H
#define KRIGINGTRANSFORMATION_H

/***************************************************************************
* idwTransformation.h
*
* Project: GPRO, v 2.0
* Purpose: Header file for class GPRO::IdwTransformation
* Author:  Wang Yujing
* E-mail:  wangyujing@lreis.ac.cn
****************************************************************************
* Copyright (c) 2019. Wang Yujing
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/
#include <vector>
#include "basicTypes.h"
#include "cellSpace.h"
#include "application.h"
#include "computeLayer.h"
#include "transformation.h"
#include "krigingOperator.h"
#include <iostream>

using namespace GPRO;

class KrigingTransformation: public Transformation<double>
{
public:
    KrigingTransformation(ComputeLayer<double> *comptLayer,KrigingOperator *krigingOperator):
    _krigingOperator(krigingOperator),Transformation(comptLayer) {}
    bool Operator(const CellCoord &coord);
private:
    KrigingOperator *_krigingOperator;
};

#endif