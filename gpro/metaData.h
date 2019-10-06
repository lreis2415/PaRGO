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
        int row;
        int column;

        int myrank;
        int processor_number;

		GDALDataType dataType;
        DomDcmpType _domDcmpType;
        SpaceDims _glbDims;
        CoordBR _MBR; /* this processor's MBR is in global coordinates */
        SpaceDims _localdims;
        //CoordBR _subMBR;
        CoordBR _localworkBR; /* this processor's workBR is in local coordinates */

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