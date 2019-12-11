#ifndef IDWTRANSFORMATION_H
#define IDWTRANSFORMATION_H

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
#include "computLayer.h"
#include "transformation.h"
#include "idwOperator.h"
#include <iostream>

using namespace GPRO;

class IdwTransformation: public Transformation<double>
{
public:
    IdwTransformation(ComputLayer<double> *comptLayer,IDWOperator *idwOperator):
    _idwOperrator(idwOperator),Transformation(comptLayer) {}
    bool Operator(const CellCoord &coord);
private:
    IDWOperator *_idwOperrator;
};

#endif