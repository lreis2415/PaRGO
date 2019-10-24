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

namespace GPRO 
{
	class MetaData 
	{
	public:
		int row;
		int column;

		int cellSize;
		string format;
		string projection;
		double* pTransform;
		int noData;

		

		int myrank;
		int processor_number;

		GDALDataType dataType;
		DomDcmpType _domDcmpType;
        SpaceDims _glbDims;
		CoordBR _MBR; /* MBR is in global coordinates */
		SpaceDims _localdims;
		//CoordBR _subMBR;
		CoordBR _localworkBR; /* workBR is in local coordinates */
	public:
		MetaData(){pTransform = new double[6];}
		~MetaData(){delete []pTransform;}
		int LoctoGloRow(int i){ return _MBR.minIRow() + i;}
		int LoctoGloCol(int j){ return _MBR.minICol() + j;}

	};
};


#endif