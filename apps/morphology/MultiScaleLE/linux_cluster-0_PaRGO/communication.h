#ifndef COMMUNICATION_H
#define COMMUNICATION_H

/***************************************************************************
* communication.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::Communication
* Author:  Zhan Lijun
* E-mail:  zhanlj@lreis.ac.cn
****************************************************************************
* Copyright (c) 2013. Zhan Lijun
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include "basicTypes.h"
#include "basicCell.h"
#include "metaData.h"
#include "neighborhood.h"
#include "rasterLayer.h"
#include <string>
#include <iostream>
#include <typeinfo>
#include <vector>
#include "mpi.h"
#include <gdal_priv.h>


namespace GPRO
{
  template <class elemType>
  class Communication
  {
	public:
		Communication();
		Communication(vector<RasterLayer<elemType>* >* pComVec);
		bool setCommunication(vector<RasterLayer<elemType>* >* pComVec);
		bool rowComm();
	private:
		vector<RasterLayer<elemType>* >* _pComVec;
		Neighborhood<elemType> * _pNbrhood;
		MetaData * _pMetadata;
  };

};

template <class elemType>
inline GPRO::Communication<elemType>::
Communication()
  :
   _pNbrhood(NULL),
   _pMetadata(NULL),
   _pComVec(NULL)
   {}

template <class elemType>
GPRO::Communication<elemType>::
Communication(vector<RasterLayer<elemType>* >* pComVec)
{
	this->_pComVec = pComVec;
	cout<<pComVec->size()<<endl;
	cout<<_pComVec->size()<<endl;
	cout<<"this->_pComVec = pComVec"<<endl;
	
	if(this->_pComVec->size() > 0)
	{
		this->_pComVec->at(0);
		cout<<"this->_pComVec->at(0);"<<endl;
		_pNbrhood = this->_pComVec->at(0)->nbrhood();
		cout<<"_pNbrhood = this->_pComVec->at(0)->nbrhood();"<<endl;
		_pMetadata =  this->_pComVec->at(0)->_pMetaData;
		cout<<"_pMetadata =  this->_pComVec->at(0)->_pMetaData;"<<endl;
	}
	cout<<"ok!"<<endl;
	
}

template <class elemType>
bool GPRO::Communication<elemType>::
setCommunication(vector<RasterLayer<elemType>* >* pComVec)
{
	_pComVec = pComVec;
	if(_pComVec->size() > 0)
	{
		_pNbrhood = _pComVec->at(0)->nbrhood();
		_pMetadata =  _pComVec->at(0)->_pMetaData;
	}
	return true;
}

template <class elemType>
bool GPRO::Communication<elemType>::
rowComm()
{
	MPI_Status status;
	for(int i = 0; i < _pComVec->size(); i++)
	{
		//cout<<"_pComVec->size is "<<_pComVec->size()<<endl;
		CellSpace<elemType>* pCellSpace = _pComVec->at(i)->cellSpace();
		MPI_Datatype datatype = _pComVec->at(i)->getMPIType();
		elemType *_matrix = pCellSpace->_matrix;
		if(_pMetadata->processor_number > 1)
		{
			int upSize = _pNbrhood->maxIRow() * _pMetadata->_glbDims.nCols();
			int lowSize = (0 - _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();

			int sendUpPoint = _pMetadata->_localworkBR.minIRow() * _pMetadata->_glbDims.nCols();
			int sendLowPoint = (_pMetadata->_localworkBR.maxIRow() + 1 + _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();

			int recvUpPoint = 0;
			int recvLowPoint = (_pMetadata->_localworkBR.maxIRow() + 1) * _pMetadata->_glbDims.nCols();

			if(_pMetadata->myrank == 0)
			{
				MPI_Send(_matrix + sendLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD);
				MPI_Recv(_matrix + recvLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD, &status);
			}
			else if( _pMetadata->myrank == ( _pMetadata->processor_number - 1) )
			{
				MPI_Send(_matrix + sendUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD);
				MPI_Recv(_matrix + recvUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD, &status);
			}
			else
			{
				MPI_Send(_matrix + sendUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD);
				MPI_Send(_matrix + sendLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD);
				MPI_Recv(_matrix + recvUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(_matrix + recvLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD, &status);

			}
		}
		
	}
	return true;
}

#endif