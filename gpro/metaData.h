/**
 * \file basicCell
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::BasicCell
 * \version 1.0
 * 
 * \copyright Copyright (c) 2013-2020
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2020 - Yujing Wang - Code reformat
 */

#ifndef METADATA_H
#define METADATA_H

#include <utility>
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

    /**
     * \ingroup gpro
     * \class MetaData 
     * \brief information to describe a raster layer.
     *  
     * In parallel situations, an actual BR is the global indices, and a virtual BR is local indices.
     *  
     * e.g. 4 processes share 80 lines.
     * 
     * No. | ActualBR | VirtualBR
     * --|--|-- 
     * 1 | 0-19 | 0-19
     * 2 | 20-39 | 0-19
     * 3 | 40-59 | 0-19
     * 4 | 60-79 | 0-19
     */
    class MetaData {
    public:
        string format; ///< format of data (GTiff etc.)
        string projection; ///< map projection
        double* pTransform; ///< map transform

        double noData; ///< noData value
        float cellSize; ///< resolution
        int row; ///< row num of data
        int column; ///< column num of data

        int myrank; ///< rank of this process
        int processor_number; ///< number of processors

        GDALDataType dataType; ///< data type of GDAL
        DomDcmpType _domDcmpType; ///< domain decomposition type. i.e. ROW, COL, BLOCK, NON
        SpaceDims _glbDims; ///< space dimension globally.
        CoordBR _MBR; ///< ACTUAL BR. this processor's MBR in global coordinates */
        SpaceDims _localdims; ///< space dimension locally
        CoordBR _localworkBR; ///< VIRTUAL BR. this processor's workBR in local coordinates

    public:
        MetaData() { pTransform = new double[6]; }

        ~MetaData() {
            if (pTransform) {
                delete pTransform;
                pTransform = nullptr;
            }
        }

        int LoctoGloRow(int i) { return _MBR.minIRow() + i; }

        int LoctoGloCol(int j) { return _MBR.minICol() + j; }

    };
};


#endif
