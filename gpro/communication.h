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
 *  - 1. 2020 - Wang Yujing - Code reformat
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
     * \brief interprocess communication
     */
    template <class elemType>
    class Communication {
    public:
        ///< Constructor
        Communication();
        ///< Construct by Communication Vector
        Communication(vector<RasterLayer<elemType>*>* pComVec);
        ///< Destructor
        ~Communication();

        ///< set Communication Vector
        bool setCommunication(vector<RasterLayer<elemType>*>* pComVec);
        ///< set MPI buffer
        bool setBuffer();
        ///< row-wise communication
        bool rowComm();
        ///< col-wise communication
        bool colComm();
    private:
        vector<RasterLayer<elemType>*>* _pComVec; ///< communication vector
        Neighborhood<elemType>* _pNbrhood; ///< neighborhood init by nbr file
        MetaData* _pMetadata;
        void* buf;
        int bufsize;
    };

};

template <class elemType>
GPRO::Communication<elemType>::
Communication()
    :
    _pComVec(nullptr),
    _pNbrhood(nullptr),
    _pMetadata(nullptr),
    buf(nullptr),
    bufsize(10) {
}

template <class elemType>
GPRO::Communication<elemType>::
Communication(vector<RasterLayer<elemType>*>* pComVec)
    :
    buf(nullptr),
    bufsize(10) {
    this->_pComVec = pComVec;
    if (this->_pComVec->size() > 0) {
        this->_pComVec->at(0);
        _pNbrhood = this->_pComVec->at(0)->nbrhood();
        _pMetadata = this->_pComVec->at(0)->_pMetaData;
    }
}


template <class elemType>
GPRO::Communication<elemType>::
~Communication() {
    if (buf != nullptr) {
        free(buf);
    }
}

template <class elemType>
bool GPRO::Communication<elemType>::
setBuffer() {
    // Buffer may be set after considering num of processes? What if there is only 1 process?
    if (_pMetadata->_domDcmpType == ROWWISE_DCMP) {
        int styleSize;
        MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize); //why 'double'? Buffer is related to the maximum data type?
        int up = _pNbrhood->maxIRow() * _pMetadata->_glbDims.nCols();
        int low = (0 - _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();
        bufsize = _pComVec->size() * (MPI_BSEND_OVERHEAD + styleSize * (up + up + low + low));
        buf = malloc(bufsize);
        MPI_Buffer_attach(buf, bufsize);
    }
    else {
        if (_pMetadata->_domDcmpType == COLWISE_DCMP) {
            int styleSize;
            MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize);
            int left = _pNbrhood->maxICol() * _pMetadata->_glbDims.nRows();
            int right = (0 - _pNbrhood->minICol()) * _pMetadata->_glbDims.nRows();
            bufsize = _pComVec->size() * (MPI_BSEND_OVERHEAD + styleSize * (left + left + right + right));
            buf = malloc(bufsize);
            MPI_Buffer_attach(buf, bufsize);
        }
        else {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: unable to support this _domDcmpType"
                << endl;
            return false;
        }
    }
    return true;
}


template <class elemType>
bool GPRO::Communication<elemType>::
setCommunication(vector<RasterLayer<elemType>*>* pComVec) {
    _pComVec = pComVec;
    if (_pComVec->size() > 0) {
        _pNbrhood = _pComVec->at(0)->nbrhood();
        _pMetadata = _pComVec->at(0)->_pMetaData;
    }
    return true;
}


template <class elemType>
bool GPRO::Communication<elemType>::
rowComm() {
    // For the moment, different layers process and communicate separately. Integrate layers into one array, and communicate once will be faster?
    MPI_Status status;
    for (std::size_t i = 0; i < _pComVec->size(); i++) {
        CellSpace<elemType>* pCellSpace = _pComVec->at(i)->cellSpace();
        MPI_Datatype datatype = _pComVec->at(i)->getMPIType();
        elemType* _matrix = pCellSpace->_matrix;
        if (_pMetadata->processor_number > 1) {
            int upSize = _pNbrhood->maxIRow() * _pMetadata->_glbDims.nCols();
            int lowSize = (0 - _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();

            int sendUpPoint = _pMetadata->_localworkBR.minIRow() * _pMetadata->_glbDims.nCols();
            int sendLowPoint = (_pMetadata->_localworkBR.maxIRow() + 1 + _pNbrhood->minIRow()) * _pMetadata->_glbDims.nCols();

            int recvUpPoint = 0;
            int recvLowPoint = (_pMetadata->_localworkBR.maxIRow() + 1) * _pMetadata->_glbDims.nCols();

            //For 3*3 neighbor, upsize==lowsize. But it may be fault here, need a asymmetric neighbor test.
            if (_pMetadata->myrank == 0) {
                MPI_Bsend(_matrix + sendLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD);
                MPI_Recv(_matrix + recvLowPoint, lowSize, datatype, 1, 1, MPI_COMM_WORLD, &status);
            }
            else if (_pMetadata->myrank == (_pMetadata->processor_number - 1)) {
                MPI_Bsend(_matrix + sendUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD);
                MPI_Recv(_matrix + recvUpPoint, upSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD, &status);
            }
            else {
                MPI_Bsend(_matrix + sendUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD);
                MPI_Bsend(_matrix + sendLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD);
                MPI_Recv(_matrix + recvUpPoint, upSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(_matrix + recvLowPoint, lowSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD, &status);
            }
        }
    }
    MPI_Buffer_detach(&buf, &bufsize);
    MPI_Buffer_attach(buf, bufsize); //why attach again?

    return true;
}


template <class elemType>
bool GPRO::Communication<elemType>::
colComm() {
    MPI_Status status;
    for (std::size_t i = 0; i < _pComVec->size(); i++) {
        CellSpace<elemType>* pCellSpace = _pComVec->at(i)->cellSpace();
        MPI_Datatype datatype = _pComVec->at(i)->getMPIType();
        elemType* _matrix = pCellSpace->_matrix;
        if (_pMetadata->processor_number > 1) {
            int leftSize = _pNbrhood->maxICol() * _pMetadata->_glbDims.nRows();
            int rightSize = (0 - _pNbrhood->minICol()) * _pMetadata->_glbDims.nRows(); // size to send to the right side.

            //Integrate discrete cells to an array
            elemType* sendLeftMatrix;
            elemType* sendRightMatrix; //data to send to the right side.
            elemType* recvLeftMatrix;
            elemType* recvRightMatrix; //data received from the right side.

            if (_pMetadata->myrank == 0) {
                sendRightMatrix = new elemType[rightSize];
                recvRightMatrix = new elemType[leftSize];
                int rightCell = 0;
                for (int cIRow = _pMetadata->_MBR.minIRow(); cIRow <= _pMetadata->_MBR.maxIRow(); ++cIRow) {
                    for (int cICol = _pMetadata->_localworkBR.maxICol() + _pNbrhood->minICol() + 1;
                         cICol <= _pMetadata->_localworkBR.maxICol(); ++cICol) {
                        sendRightMatrix[rightCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++rightCell;
                    }
                }
            }
            else if (_pMetadata->myrank == (_pMetadata->processor_number - 1)) {
                sendLeftMatrix = new elemType[leftSize];
                recvLeftMatrix = new elemType[rightSize];
                int leftCell = 0;
                for (int cIRow = _pMetadata->_MBR.minIRow(); cIRow <= _pMetadata->_MBR.maxIRow(); ++cIRow) {
                    for (int cICol = _pMetadata->_localworkBR.minICol();
                         cICol <= _pMetadata->_localworkBR.minICol() + _pNbrhood->maxICol() - 1; ++cICol) {
                        sendLeftMatrix[leftCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++leftCell;
                    }
                }
            }
            else {
                sendLeftMatrix = new elemType[leftSize];
                sendRightMatrix = new elemType[rightSize];
                recvLeftMatrix = new elemType[rightSize];
                recvRightMatrix = new elemType[leftSize];
                int leftCell = 0, rightCell = 0;
                for (int cIRow = _pMetadata->_MBR.minIRow(); cIRow <= _pMetadata->_MBR.maxIRow(); ++cIRow) {
                    for (int cICol = _pMetadata->_localworkBR.minICol();
                         cICol <= _pMetadata->_localworkBR.minICol() + _pNbrhood->maxICol() - 1; ++cICol) {
                        sendLeftMatrix[leftCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++leftCell;
                    }
                    for (int cICol = _pMetadata->_localworkBR.maxICol() + _pNbrhood->minICol() + 1;
                         cICol <= _pMetadata->_localworkBR.maxICol(); ++cICol) {
                        sendRightMatrix[rightCell] = _matrix[cIRow * _pMetadata->_localdims.nCols() + cICol];
                        ++rightCell;
                    }
                }
            }

            //communicate
            if (_pMetadata->myrank == 0) {
                MPI_Bsend(sendRightMatrix, rightSize, datatype, 1, 1, MPI_COMM_WORLD);
                MPI_Recv(recvRightMatrix, leftSize, datatype, 1, 1, MPI_COMM_WORLD, &status);
            }
            else if (_pMetadata->myrank == (_pMetadata->processor_number - 1)) {
                MPI_Bsend(sendLeftMatrix, leftSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD);
                MPI_Recv(recvLeftMatrix, rightSize, datatype, _pMetadata->processor_number - 2, 1, MPI_COMM_WORLD, &status);
            }
            else {
                MPI_Bsend(sendLeftMatrix, leftSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD);
                MPI_Bsend(sendRightMatrix, rightSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD);
                MPI_Recv(recvLeftMatrix, rightSize, datatype, _pMetadata->myrank - 1, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(recvRightMatrix, leftSize, datatype, _pMetadata->myrank + 1, 1, MPI_COMM_WORLD, &status);
            }

            //update cells according to receive arrays
            //Is method in the CellSpace class able to do the following work?
            int cRow = 0, cCol = 0;
            if (_pMetadata->myrank == 0) {
                for (int i = 0; i < leftSize; ++i) {
                    cRow = i / _pNbrhood->maxICol() + _pMetadata->_MBR.minIRow();
                    cCol = i % _pNbrhood->maxICol() + _pMetadata->_localworkBR.maxICol() + 1;
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvRightMatrix[i];
                }
                delete[]sendRightMatrix; //why it doesn't allow delete here, error occurs
                delete[]recvRightMatrix;
            }
            else if (_pMetadata->myrank == (_pMetadata->processor_number - 1)) {
                for (int i = 0; i < rightSize; ++i) {
                    cRow = i / (0 - _pNbrhood->minICol()) + _pMetadata->_MBR.minIRow();
                    cCol = i % (0 - _pNbrhood->minICol()) + _pMetadata->_localworkBR.minICol() + _pNbrhood->minICol();
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvLeftMatrix[i];
                }
                delete[]sendLeftMatrix;
                delete[]recvLeftMatrix;
            }
            else {
                for (int i = 0; i < leftSize; ++i) {
                    cRow = i / _pNbrhood->maxICol() + _pMetadata->_MBR.minIRow();
                    cCol = i % _pNbrhood->maxICol() + _pMetadata->_localworkBR.maxICol() + 1;
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvRightMatrix[i];
                }
                for (int i = 0; i < rightSize; ++i) {
                    cRow = i / (0 - _pNbrhood->minICol()) + _pMetadata->_MBR.minIRow();
                    cCol = i % (0 - _pNbrhood->minICol()) + _pMetadata->_localworkBR.minICol() + _pNbrhood->minICol();
                    _matrix[cRow * _pMetadata->_localdims.nCols() + cCol] = recvLeftMatrix[i];
                }

                delete[]sendRightMatrix;
                delete[]recvRightMatrix;
                delete[]sendLeftMatrix;
                delete[]recvLeftMatrix;
            }
        }
    }

    MPI_Buffer_detach(&buf, &bufsize);
    MPI_Buffer_attach(buf, bufsize);
    return true;
}

#endif
