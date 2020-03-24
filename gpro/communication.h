/**
 * \file basicCell
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::Communication
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

#ifndef COMMUNICATION_H
#define COMMUNICATION_H

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

namespace GPRO {
    /**
     * \ingroup gpro
     * \class Communication
     * \brief 
     */
    template<class elemType>
    class Communication {
    public:
        /// Constructor
        Communication();
        /// Construct by Communication Vector
        Communication( vector<RasterLayer < elemType> * > *pComVec);
        /// Deconstructor
        ~Communication();

        /// set Communication Vector
        bool setCommunication( vector<RasterLayer < elemType> * > *pComVec);
        /// set MPI buffer
        bool setBuffer();
        /// row-wise communication
        bool rowComm();
        /// col-wise communication
        bool colComm();
    private:
        vector<RasterLayer < elemType>* > *_pComVec; /// communication vector
        Neighborhood <elemType> *_pNbrhood; /// neighborhood init by nbr file
        MetaData *_pMetadata;
        void *buf;
        int bufsize;
    };

};

template<class elemType>
inline GPRO::Communication<elemType>::
Communication()
    :
    _pNbrhood( NULL ),
    _pMetadata( NULL ),
    _pComVec( NULL ),
    buf( NULL ),
    bufsize( 10 ) {}

template<class elemType>
GPRO::Communication<elemType>::
Communication( vector<RasterLayer < elemType> * > *pComVec)
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
		this->_pComVec->at(0);
		//cout<<"this->_pComVec->at(0);"<<endl;
		_pNbrhood = this->_pComVec->at( 0 )->nbrhood();
		//cout<<"_pNbrhood = this->_pComVec->at(0)->nbrhood();"<<endl;
		_pMetadata = this->_pComVec->at( 0 )->_pMetaData;
		//cout<<"_pMetadata =  this->_pComVec->at(0)->_pMetaData;"<<endl;
	}
}


template<class elemType>
inline GPRO::Communication<elemType>::
~Communication() {
    //cout<<"buf is "<<buf<<endl;
    if ( buf != NULL ) {
        free( buf );
    }
}

template<class elemType>
bool GPRO::Communication<elemType>::
setBuffer() {
	/*是不是需要考虑总进程数再设置这里的buffer，如果是单进程怎么处理*/
    if ( _pMetadata->_domDcmpType == ROWWISE_DCMP ) {
        int styleSize;
        MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize );    //为什么这里要用double，什么含义；根据最大数据类型申请缓冲区大小？
        int up = _pNbrhood->maxIRow() * _pMetadata->_glbDims.nCols();
        int low = ( 0 - _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();
        bufsize = _pComVec->size() * ( MPI_BSEND_OVERHEAD + styleSize * ( up + up + low + low ));
        buf = malloc( bufsize );
        MPI_Buffer_attach( buf, bufsize );
    } else {
		if ( _pMetadata->_domDcmpType == COLWISE_DCMP ) {
			int styleSize;
			MPI_Pack_size( 1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize );
			int left = _pNbrhood->maxICol() * _pMetadata->_glbDims.nRows();
			int right = ( 0 - _pNbrhood->minICol()) * _pMetadata->_glbDims.nRows();
			bufsize = _pComVec->size() * ( MPI_BSEND_OVERHEAD + styleSize * ( left + left + right + right ));
			buf = malloc( bufsize );
			MPI_Buffer_attach( buf, bufsize );
		} else {
			cerr << __FILE__ << " " << __FUNCTION__ \
				<< " Error: unable to support this _domDcmpType" \
				<< endl;
            return false;
		}
	}
    return true;
}


template<class elemType>
bool GPRO::Communication<elemType>::
setCommunication( vector<RasterLayer < elemType> * > *pComVec)
{
	_pComVec = pComVec;
	if(_pComVec->size() > 0)
	{
		_pNbrhood = _pComVec->at( 0 )->nbrhood();
		_pMetadata = _pComVec->at( 0 )->_pMetaData;
	}
	return true;
}


template<class elemType>
bool GPRO::Communication<elemType>::
rowComm() {
	//现在是多个图层分别整合，分别通信；可以考虑多个图层整合成一个数组，一次性通信,不知能否提升效率；其他划分方式也一样
    MPI_Status status;
    for ( std::size_t i = 0; i < _pComVec->size(); i++ ) {
        //cout<<"_pComVec->size is "<<_pComVec->size()<<endl;
        CellSpace<elemType> *pCellSpace = _pComVec->at( i )->cellSpace();
        MPI_Datatype datatype = _pComVec->at( i )->getMPIType();
        elemType *_matrix = pCellSpace->_matrix;
        if ( _pMetadata->processor_number > 1 ) {
            int upSize = _pNbrhood->maxIRow() * _pMetadata->_glbDims.nCols();
            int lowSize = ( 0 - _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();

            int sendUpPoint = _pMetadata->_localworkBR.minIRow() * _pMetadata->_glbDims.nCols();
            int sendLowPoint = ( _pMetadata->_localworkBR.maxIRow() + 1 + _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();

            int recvUpPoint = 0;
            int recvLowPoint = ( _pMetadata->_localworkBR.maxIRow() + 1 ) * _pMetadata->_glbDims.nCols();

            /*对目前的3*3邻域来说，这里的upsize==lowsize，但是怀疑这里写错了，待写个不对称的邻域测试*/
            if ( _pMetadata->myrank == 0 ) {
                MPI_Bsend( _matrix + sendLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD );
                MPI_Recv( _matrix + recvLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD, &status );
            } else if ( _pMetadata->myrank == ( _pMetadata->processor_number - 1 )) {
                MPI_Bsend( _matrix + sendUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD );
                MPI_Recv( _matrix + recvUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD, &status );
            } else {
                MPI_Bsend( _matrix + sendUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD );
                MPI_Bsend( _matrix + sendLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD );
                MPI_Recv( _matrix + recvUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD, &status );
                MPI_Recv( _matrix + recvLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD, &status );
            }
        }
    }
    MPI_Buffer_detach( &buf, &bufsize );
    MPI_Buffer_attach( buf, bufsize );	//为什么还要再申请，不应该有吧？

    return true;
}


template<class elemType>
bool GPRO::Communication<elemType>::
colComm() {
    MPI_Status status;
    for ( std::size_t i = 0; i < _pComVec->size(); i++ ) {
        CellSpace<elemType> *pCellSpace = _pComVec->at( i )->cellSpace();
        MPI_Datatype datatype = _pComVec->at( i )->getMPIType();
        elemType *_matrix = pCellSpace->_matrix;
        if ( _pMetadata->processor_number > 1 ) {
            int leftSize = _pNbrhood->maxICol() * _pMetadata->_glbDims.nRows();
            int rightSize = ( 0 - _pNbrhood->minICol()) * _pMetadata->_glbDims.nRows();//向右发送的大小

            //Integrate discrete cells to an array
            elemType *sendLeftMatrix;
            elemType *sendRightMatrix;    //向右发送的数据
            elemType *recvLeftMatrix;
            elemType *recvRightMatrix;    //从右边接收来的数据

            if ( _pMetadata->myrank == 0 ) {
                sendRightMatrix = new elemType[rightSize];
                recvRightMatrix = new elemType[leftSize];
                int rightCell = 0;
                for ( int cIRow = _pMetadata->_MBR.minIRow(); cIRow <= _pMetadata->_MBR.maxIRow(); ++cIRow ) {
                    for ( int cICol = _pMetadata->_localworkBR.maxICol() + _pNbrhood->minICol() + 1;
                          cICol <= _pMetadata->_localworkBR.maxICol(); ++cICol ) {
                        sendRightMatrix[rightCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++rightCell;
                    }
                }
            } else if ( _pMetadata->myrank == ( _pMetadata->processor_number - 1 )) {
                sendLeftMatrix = new elemType[leftSize];
                recvLeftMatrix = new elemType[rightSize];
                int leftCell = 0;
                for ( int cIRow = _pMetadata->_MBR.minIRow(); cIRow <= _pMetadata->_MBR.maxIRow(); ++cIRow ) {
                    for ( int cICol = _pMetadata->_localworkBR.minICol();
                          cICol <= _pMetadata->_localworkBR.minICol() + _pNbrhood->maxICol() - 1; ++cICol ) {
                        sendLeftMatrix[leftCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++leftCell;
                    }
                }
            } else {
                sendLeftMatrix = new elemType[leftSize];
                sendRightMatrix = new elemType[rightSize];
                recvLeftMatrix = new elemType[rightSize];
                recvRightMatrix = new elemType[leftSize];
                int leftCell = 0, rightCell = 0;
                for ( int cIRow = _pMetadata->_MBR.minIRow(); cIRow <= _pMetadata->_MBR.maxIRow(); ++cIRow ) {
                    for ( int cICol = _pMetadata->_localworkBR.minICol();
                          cICol <= _pMetadata->_localworkBR.minICol() + _pNbrhood->maxICol() - 1; ++cICol ) {
                        sendLeftMatrix[leftCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++leftCell;
                    }
                    for ( int cICol = _pMetadata->_localworkBR.maxICol() + _pNbrhood->minICol() + 1;
                          cICol <= _pMetadata->_localworkBR.maxICol(); ++cICol ) {
                        sendRightMatrix[rightCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++rightCell;
                    }
                }
            }

            //communicate
            if ( _pMetadata->myrank == 0 ) {
                MPI_Bsend( sendRightMatrix, rightSize, datatype, 1, 1, MPI_COMM_WORLD );
                MPI_Recv( recvRightMatrix, leftSize, datatype, 1, 1, MPI_COMM_WORLD, &status );
            } else if ( _pMetadata->myrank == ( _pMetadata->processor_number - 1 )) {
                MPI_Bsend( sendLeftMatrix, leftSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD );
                MPI_Recv( recvLeftMatrix, rightSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD, &status );
            } else {
                MPI_Bsend( sendLeftMatrix, leftSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD );
                MPI_Bsend( sendRightMatrix, rightSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD );
                MPI_Recv( recvLeftMatrix, rightSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD, &status );
                MPI_Recv( recvRightMatrix, leftSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD, &status );
            }
            
			//update cells according to receive arrays
			/*这一步是不是利用cellcpace类中封装的功能就可实现*/
            int cRow = 0, cCol = 0;
            if ( _pMetadata->myrank == 0 ) {
                for ( int i = 0; i < leftSize; ++i ) {
                    cRow = i / _pNbrhood->maxICol() + _pMetadata->_MBR.minIRow();
                    cCol = i % _pNbrhood->maxICol() + _pMetadata->_localworkBR.maxICol() + 1;
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvRightMatrix[i];
                }
                delete[]sendRightMatrix;    //为什么这里不允许释放内存,会出错
                delete[]recvRightMatrix;
            } else if ( _pMetadata->myrank == ( _pMetadata->processor_number - 1 )) {
                for ( int i = 0; i < rightSize; ++i ) {
                    cRow = i / ( 0 - _pNbrhood->minICol()) + _pMetadata->_MBR.minIRow();
                    cCol = i % ( 0 - _pNbrhood->minICol()) + _pMetadata->_localworkBR.minICol() + _pNbrhood->minICol();
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvLeftMatrix[i];
                }
                delete[]sendLeftMatrix;
                delete[]recvLeftMatrix;
            } else {
                for ( int i = 0; i < leftSize; ++i ) {
                    cRow = i / _pNbrhood->maxICol() + _pMetadata->_MBR.minIRow();
                    cCol = i % _pNbrhood->maxICol() + _pMetadata->_localworkBR.maxICol() + 1;
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvRightMatrix[i];
                }
                for ( int i = 0; i < rightSize; ++i ) {
                    cRow = i / ( 0 - _pNbrhood->minICol()) + _pMetadata->_MBR.minIRow();
                    cCol = i % ( 0 - _pNbrhood->minICol()) + _pMetadata->_localworkBR.minICol() + _pNbrhood->minICol();
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvLeftMatrix[i];
                }

                delete[]sendRightMatrix;
                delete[]recvRightMatrix;
                delete[]sendLeftMatrix;
                delete[]recvLeftMatrix;
            }
        }
    }

    MPI_Buffer_detach( &buf, &bufsize );
    MPI_Buffer_attach( buf, bufsize );
    //cout << _pMetadata->myrank << " done." << endl;
    return true;
}

#endif