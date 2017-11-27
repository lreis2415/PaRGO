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
		~Communication();
		bool setCommunication(vector<RasterLayer<elemType>* >* pComVec);
		bool rowComm();
		bool setBuffer();
	private:
		vector<RasterLayer<elemType>* >* _pComVec;
		Neighborhood<elemType> * _pNbrhood;
		MetaData * _pMetadata;
		void* buf;
		int bufsize;
  };

};

template <class elemType>
inline GPRO::Communication<elemType>::
Communication()
  :
   _pNbrhood(NULL),
   _pMetadata(NULL),
   _pComVec(NULL),
   buf(NULL),
   bufsize(10)
   {}

template <class elemType>
GPRO::Communication<elemType>::
Communication(vector<RasterLayer<elemType>* >* pComVec)
	:
	buf(NULL),
	bufsize(10)
{
	this->_pComVec = pComVec;
	//cout<<pComVec->size()<<endl;
	//cout<<_pComVec->size()<<endl;
	//cout<<"this->_pComVec = pComVec"<<endl;
	
	if(this->_pComVec->size() > 0)
	{
		//this->_pComVec->at(0);
		//cout<<"this->_pComVec->at(0);"<<endl;
		_pNbrhood = this->_pComVec->at(0)->nbrhood();
		//cout<<"_pNbrhood = this->_pComVec->at(0)->nbrhood();"<<endl;
		_pMetadata =  this->_pComVec->at(0)->_pMetaData;
		//cout<<"_pMetadata =  this->_pComVec->at(0)->_pMetaData;"<<endl;
	}
	//cout<<"ok!"<<endl;
	
}

template <class elemType>
inline GPRO::Communication<elemType>::
~Communication()
{
	if(buf != NULL)
	{
		free(buf);
	}
}

template <class elemType>
bool GPRO::Communication<elemType>::
setBuffer()
{
	if( _pComVec->size() >0 ){
		int styleSize;
		MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize );
		int up = _pNbrhood->maxIRow() * _pMetadata->_glbDims.nCols();
		int low = (0 - _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();
		bufsize = _pComVec->size() * ( MPI_BSEND_OVERHEAD + styleSize*(up+up+low+low) );
		buf = malloc( bufsize );
		MPI_Buffer_attach( buf, bufsize );
	}
	return true;
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
	if( _pComVec->size() == 0 ){
		return true;
	}
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
				MPI_Bsend(_matrix + sendLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD);
				MPI_Recv(_matrix + recvLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD, &status);
			}
			else if( _pMetadata->myrank == ( _pMetadata->processor_number - 1) )
			{
				MPI_Bsend(_matrix + sendUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD);
				MPI_Recv(_matrix + recvUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD, &status);
			}
			else
			{
				MPI_Bsend(_matrix + sendUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD);
				MPI_Bsend(_matrix + sendLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD);
				MPI_Recv(_matrix + recvUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(_matrix + recvLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD, &status);

			}
		}
		
	}
	MPI_Buffer_detach(&buf,&bufsize);
	MPI_Buffer_attach(buf,bufsize);
	return true;
}

#endif