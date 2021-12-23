/**
 * \file transformation
 * \author Ai Beibei (aibb@lreis.ac.cn)
 * \brief Header file for class GPRO::Transformation
 * \version 2.0
 * 
 * \copyright Copyright (c) 2013-2020
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2020 - Wang Yujing - Code reformat
 */

#ifndef Transformation_H
#define Transformation_H

#include <vector>
#include "basicTypes.h"
#include "cellSpace.h"
#include "application.h"
#include "computeLayer.h"
#include <iostream>

using namespace std;

namespace GPRO {
    /**
     * \ingroup gpro
     * \class Transformation 
     * \brief Transformation from data domain (RasterLayer) to compute domain (ComputeLayer)
     */
    template <class elemType>
    class Transformation {
    public:
        Transformation();

        /**
         * \brief To be overloaded by subclass.
         * \param[in] pLayer the compute layer to be derived.
         */
        Transformation(ComputeLayer<elemType>* pLayer);

        /**
         * \brief A default implementation for uniform distribution of data but nonuniform distribution of computation.
         * \param[in] nodataLoad intensity value for NoData grid.
         * \param[in] validLoad intensity value for valid grid.
         * \param[in] pLayer the compute layer to be derived.
         */
        Transformation(elemType nodataLoad, elemType validLoad, ComputeLayer<elemType>* pLayer);

        virtual ~Transformation() {
        } //virtual Destructor

        virtual bool isTermination(bool isIter = false) { return isIter; } // confusing.
        bool Configure(ComputeLayer<elemType>* pLayer, bool isCommunication);
        bool paramInit();
        virtual bool Operator(const CellCoord& coord);

        /**
         * \brief To start the process.
         */
        bool run();

    protected:
        elemType _nodataLoad; ///< intensity value for NoData grid
        elemType _validLoad; ///< intensity value for valid grid
        int _myRank; ///< rank of this process

        double _noData; ///< empty value in data domain
        double _computGrain; ///< granularity for compute layer. 1 compute layer grid contains _computGrain^2 raster layer grids.
        CoordBR _dataMBR; ///< MBR of data domain
        CoordBR _dataWorkBR; ///< workBR of this process

        //vector<RasterLayer<elemType>* > CommVec;	// to be decided. It's needed for parallel transformation
        ComputeLayer<elemType>* _pComptLayer; ///< the output
        //vector<RasterLayer<elemType>* > _pDataLayersV;	// to implement some default versions
        CoordBR* _pWorkBR; ///< workBR of this Transformation
        int Termination; ///< if it terminates after this round of traverse
    };
};

template <class elemType>
GPRO::Transformation<elemType>::
Transformation()
    : _nodataLoad(0),
      _validLoad(0),
      _pComptLayer(nullptr),
      _pWorkBR(nullptr),
      Termination(1) {
}


template <class elemType>
GPRO::Transformation<elemType>::
Transformation(ComputeLayer<elemType>* pLayer)
    : _nodataLoad(0),
      _validLoad(0),
      _pComptLayer(pLayer),
      _pWorkBR(nullptr),
      Termination(1) {
    Configure(pLayer, false);
    if (pLayer->_pMetaData != NULL) {
        _myRank = pLayer->metaData()->myrank;
    }
}


template <class elemType>
GPRO::Transformation<elemType>::
Transformation(elemType nodataLoad, elemType validLoad, ComputeLayer<elemType>* pLayer)
    : _nodataLoad(nodataLoad),
      _validLoad(validLoad),
      _pComptLayer(pLayer),
      _pWorkBR(nullptr),
      Termination(1) {
    Configure(pLayer, false);
    if (pLayer->_pMetaData != NULL) {
        _myRank = pLayer->_pMetaData->myrank;
    }
}


template <class elemType>
bool GPRO::Transformation<elemType>::
Configure(ComputeLayer<elemType>* pLayer, bool isCommunication) {
    //It is simple for the moment, may have more members or communications yet.
    if (_pWorkBR == nullptr && &pLayer->_pMetaData != NULL) {
        _pWorkBR = &pLayer->_pMetaData->_localworkBR;
    }
    if (isCommunication) {
        //do sth.
    }

    return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
paramInit() {
    if (_pComptLayer->_vDataLayers.empty()) {
        cerr << "Datalayers used for calculating compute layer should not be null." << endl;
        return false;
    }
    _noData = _pComptLayer->_vDataLayers[0]->_pMetaData->noData;
    _computGrain = _pComptLayer->getComputeGrain();
    //_computGrain = (double)(_pComptLayer->_pMetaData->cellSize / _pComptLayer->_vDataLayers[0]->_pMetaData->cellSize);
    _dataMBR = _pComptLayer->_vDataLayers[0]->_pMetaData->_MBR;
    _dataWorkBR = _pComptLayer->_vDataLayers[0]->_pMetaData->_localworkBR;

    return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
Operator(const CellCoord& coord) {
    int cRow = coord.iRow();
    int cCol = coord.iCol();
    CellSpace<elemType>& computL = *(_pComptLayer->cellSpace());

    if (cRow == _pWorkBR->minIRow() && cCol == _pWorkBR->minICol()) {
        cout << "Transformation operator() function called." << endl;
        // if( _nodataLoad == _validLoad ){
        // 	cerr<<"The load is balanced. No need to use compute layer."<<endl;
        // 	return false;
        // }
        if (_nodataLoad < 0 || _validLoad < 0) {
            cerr << "The load specified cannot be negative." << endl;
            return false;
        }

        if (!paramInit()) {
            return false;
        }
    }
    computL[cRow][cCol] = 0.0;

    for (typename vector<RasterLayer<elemType>*>::iterator iter = _pComptLayer->_vDataLayers.begin(); iter != _pComptLayer->_vDataLayers.end(); ++iter) {
        //iterates each layer, add up to the computeLayer
        const CellSpace<elemType>& dataL = *(*iter)->cellSpace();
        for (int dRow = cRow * _computGrain + _dataWorkBR.minIRow(); dRow < (cRow + 1) * _computGrain + _dataWorkBR.minIRow(); ++dRow) {
            for (int dCol = cCol * _computGrain + _dataWorkBR.minICol(); dCol < (cCol + 1) * _computGrain + _dataWorkBR.minICol(); ++dCol) {
                if (dRow > _dataMBR.maxIRow() || dRow >= dataL.nRows() || dCol > _dataMBR.maxICol() || dCol >= dataL.nCols()) {
                    continue;
                }
                if (fabs(dataL[dRow][dCol] - (*iter)->metaData()->noData) > Eps) {
                    computL[cRow][cCol] += _validLoad;
                }
                else {
                    computL[cRow][cCol] += _nodataLoad;
                }
            }
        }
    }
    return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
run() {
    if (Application::_programType != MPI_Type && Application::_programType != MPI_OpenMP_Type) {
        cerr << "not supported yet." << endl;
    }
    bool flag = true;
    if (GetRank() == 0) {
        //It's serial for the moment.
        int termSum = 1;
        do {
            Termination = 1;
            for (int iRow = _pWorkBR->minIRow(); iRow <= _pWorkBR->maxIRow(); iRow++) {
                for (int iCol = _pWorkBR->minICol(); iCol <= _pWorkBR->maxICol(); iCol++) {
                    CellCoord coord(iRow, iCol);
                    if (!Operator(coord)) {
                        cout << "Operator is not successes!" << endl;
                        flag = false;
                        break;
                    }
                    // printf("[%i,%i]",iRow,iCol);
                }
            }
            //MPI_Allreduce(&Termination, &termSum, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
            termSum = Termination; //It's serial for the moment.
        }
        while (!termSum);
    }
    return flag;
}


#endif
