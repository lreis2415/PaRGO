/**
 * \file basicCell
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::BasicCell
 * \version 1.0
 * 
 * \copyright Copyright (c) 2013
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2019-10 - Yujing Wang - Code reformat
 */

#ifndef METADATA_H
#define METADATA_H

#include <utility>	//这些头文件都是必须的吗，测试是否可删
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <gdal_priv.h>

using namespace std;

namespace GPRO {
    class MetaData {
    public:
        string format;
        string projection;
        double *pTransform;

		//int noData;
		double noData;
        float cellSize;
        int row; /// row num of data
        int column; /// column num of data

        int myrank; /// rank of this process
        int processor_number; /// number of processors

		GDALDataType dataType; /// data type of GDAL
        DomDcmpType _domDcmpType; /// domain decomposition type
        SpaceDims _glbDims; /// space dimension globally
        CoordBR _MBR; /* this processor's MBR in global coordinates */
        SpaceDims _localdims; /// space dimension locally
        CoordBR _localworkBR; /// this processor's workBR in local coordinates

	public:
        MetaData() { pTransform = new double[6]; }

        ~MetaData() {
            if ( pTransform ) {
                delete pTransform;
                pTransform = 0;
            }
        }
        //~MetaData();

        int LoctoGloRow( int i ) { return _MBR.minIRow() + i; }

        int LoctoGloCol( int j ) { return _MBR.minICol() + j; }

    };
};


#endif