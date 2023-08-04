/**
 * \file basicTypes
 * \author Zhan Lijun (zhanlj@lreis.ac.cn), Wang Yujing
 * \brief Header file for class GPRO::RasterOperator
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
#ifndef RASTEROPERATOR_H
#define RASTEROPERATOR_H

#include "basicTypes.h"
#include "rasterLayer.h"
#include "communication.h"
#include "application.h"
#include <iostream>
#include <typeinfo>
#include <vector>
#include "mpi.h"
#include <omp.h>
#include <gdal_priv.h>

using namespace std;

namespace GPRO {
    /**
     * \ingroup gpro
     * \class RasterOperator
     * \brief A basic super class that each Operator of a specific algorithm should extends
     */
    template <class elemType>
    class RasterOperator {
    public:
        RasterOperator()
            : _domDcmpType(NON_DCMP),
              _pComptLayer(nullptr),
              _pWorkBR(nullptr),
              commFlag(false),
              Termination(true),
              computeTimeExceptLastCell(0),
              _writePreExpLoad(false) {
        }

        virtual ~RasterOperator() {
        }


        /**
         * \brief Basic function in which serial-style algorithm should be writen.
         * \param[in] coord coordinate of the central cell(point)
         * \param[in] operFlag an unused controlling variable.
         */
        virtual bool Operator(const CellCoord& coord, bool operFlag) { return operFlag; }

        /**
         * \brief Implement to control the termination.
         */
        virtual bool isTermination() { return false; }

        /**
         * \brief Operator is invoked here. Traverse throughout the workBR.
         * \param[in] pWorkBR the Bounding Rectangle area to be traversed.
         */
        bool Work(const CoordBR* pWorkBR);

        /**
         * \brief Assign the workBR in pLayer to this Operator.
         * \param[in] pWorkBR the rectangle area to be traversed.
         * \param[in] pWorkBR the rectangle area to be traversed.
         */
        bool Configure(RasterLayer<elemType>* pLayer, bool isCommunication);

        /**
         * \brief This entrance method is to be called in the main process.
         */
        bool Run();

        /**
         * \brief Mount the compute layer to this Operator.
         *  
         * It is designed for intensity evaluation strategy. The member comptLayer will be initialised.
         * \param[in] layerD the rectangle area to be traversed.
         * \param[in] pWorkBR the rectangle area to be traversed.
         */
        void comptLayer(RasterLayer<elemType>& layerD);

    private:
        DomDcmpType _domDcmpType;

    public:
        RasterLayer<elemType>* _pComptLayer; ///< compute layer, mounted to be filled
        // char* _computeLayerOutPath;

        vector<RasterLayer<elemType>*> CommVec; ///< raster layers to communicate between processes
        CoordBR* _pWorkBR; ///< work bounding rectangle of this process
        bool commFlag; ///< true if involve communication
        int Termination; ///< typically 1 implies terminate, 0 implies another traversion
        double computeTimeExceptLastCell;
        double operatorReduceTime;
        bool _writePreExpLoad;
    };
};

template <class elemType>
void GPRO::RasterOperator<elemType>::
comptLayer(RasterLayer<elemType>& layerD) {
    _pComptLayer = &layerD;
    _pComptLayer->cellSpace()->initVals(0); //wyj 2019-12-7 init for fcm idw Operator
    //Configure(_pComptLayer, false); // wyj 2020-7-29 the overwriting of workBR from different layer size (data layer greater and compute layer smaller) causes error
}


// template <class elemType>
// void GPRO::RasterOperator<elemType>::
// comptLayer(RasterLayer<elemType>& layerD,char *computeLayerOutPath) {
//     comptLayer(layerD);
//     _computeLayerOutPath=computeLayerOutPath;
// }

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Configure(RasterLayer<elemType>* pLayer, bool isCommunication) {
    //wyj: why not overwrite workBR? It's OK tentatively.
    _pWorkBR = &pLayer->_pMetaData->_localworkBR;
    _domDcmpType = pLayer->_pMetaData->_domDcmpType;

    //if (_pWorkBR == NULL) { //wyj: why not overwrite workBR?
    //    _pWorkBR = &pLayer->_pMetaData->_localworkBR;
    //    _domDcmpType = pLayer->_pMetaData->_domDcmpType;
    //}
    if (isCommunication) {
        CommVec.push_back(pLayer);
        if (commFlag == false) {
            commFlag = true;
        }
    }

    return true;
}

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Work(const CoordBR* const pWBR) {
    bool flag = true; //if this method finishes normally

    //cout<<"pWBR->minIRow() "<<pWBR->minIRow()<<"  pWBR->maxIRow()"<<pWBR->maxIRow()<<endl;
    //cout<<"_pCommVec->size() "<<CommVec.size()<<endl;
    Communication<elemType> COMNI(&CommVec);
    if (commFlag) {
        COMNI.setBuffer();
    }
    int noterm = 1; // if iteration continues. "may change to other vars." (?)
    int itera = 0; // iteration nums
    int myRank = 0;
    int nRow = pWBR->maxIRow() - pWBR->minIRow();
    int delim = 20;
    int lastRowInterval = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    double startTime = MPI_Wtime();
    double endTime = 0;
    double iterStartTime = 0;
    double intervalStartTime = MPI_Wtime();
    if (Application::_programType == MPI_Type) {
        MPI_Barrier(MPI_COMM_WORLD);
        do {
            iterStartTime = MPI_Wtime();
            Termination = 1;
            for (int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) {
                for (int iCol = pWBR->minICol(); iCol <= pWBR->maxICol(); iCol++) {
                    CellCoord coord(iRow, iCol);
                    //if(iRow==pWBR->maxIRow()&&iCol== pWBR->maxICol()) {
                    //    computeTimeExceptLastCell+=MPI_Wtime()-iterStartTime;
                    //}
                    if (!Operator(coord, true)) {
                        cout << "Operator is not successes!" << endl;
                        flag = false;
                        break;
                    }
                }
#ifdef _DEBUG
                int rowInterval = (iRow - pWBR->minIRow()) / (static_cast<double>(nRow) / delim);
                if (rowInterval != lastRowInterval) {
                    lastRowInterval = rowInterval;
                    cout << "rank" << myRank << "\t[";
                    for (int i = 0; i < lastRowInterval; ++i) {
                        cout << ".";
                    }
                    for (int i = 0; i < delim - lastRowInterval; ++i) {
                        cout << " ";
                    }
                    endTime = MPI_Wtime();
                    cout << "]" << endTime - intervalStartTime << "s (" << iRow - nRow / delim << "~" << iRow << ")" << endl;
                    intervalStartTime = MPI_Wtime();
                }
#endif
            }
            //cout<<"rank"<<myRank<<" compute time "<<computeTimeExceptLastCell<<"s"<<endl;
            //cout<<"rank"<<myRank<<" reduce time "<<operatorReduceTime<<"s"<<endl;
            //cout<<"rank"<<myRank<<" iter time "<<MPI_Wtime()-iterStartTime<<"s"<<endl;
            if (commFlag) {
                if (_domDcmpType == ROWWISE_DCMP) {
                    COMNI.rowComm();
                }
                if (_domDcmpType == COLWISE_DCMP) {
                    COMNI.colComm();
                }
            }
            MPI_Allreduce(&Termination, &noterm, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
            itera++;
        }
        while (!noterm);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else if (Application::_programType == MPI_OpenMP_Type) {
        MPI_Barrier(MPI_COMM_WORLD);
        do {
            Termination = 1;
#pragma omp parallel for
            for (int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) {
                for (int iCol = pWBR->minICol(); iCol <= pWBR->maxICol(); iCol++) {
                    CellCoord coord(iRow, iCol);
                    if (!Operator(coord, true)) {
                        cout << "Operator is not sucessess!" << endl;
                        flag = false;
                        break;
                    }
                }
            }

            if (commFlag) {
                if (_domDcmpType == ROWWISE_DCMP) {
                    COMNI.rowComm();
                }
                if (_domDcmpType == COLWISE_DCMP) {
                    COMNI.colComm();
                }
            }

            MPI_Allreduce(&Termination, &noterm, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

            itera++;

        }
        while (!noterm);
        MPI_Barrier(MPI_COMM_WORLD);
        cout << "iterative numerber is " << itera << endl;
        //		MPI_Barrier(MPI_COMM_WORLD);
        //		double tcount=0;
        //		double tcount2=0;
        //		do
        //		{
        //			Termination = 1;
        //			double starttime;
        //			double endtime;
        //			double starttime2;
        //			double endtime2;
        //
        //			if (itera==0)
        //			{
        //				processing();
        //			}
        //			starttime = MPI_Wtime();
        //
        //#pragma omp parallel 
        //			{
        //			//#pragma omp for schedule(dynamic,8)
        //#pragma omp for
        //				for(int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) 
        //				{
        //					for(int iCol = pWBR->minICol();	iCol <= pWBR->maxICol(); iCol++) 
        //					{
        //						CellCoord coord(iRow, iCol);
        //						if(!Operator(coord)) 
        //						{
        //							cout<<"Operator is not sucessess!"<<endl;
        //							flag = false;
        //							break;
        //
        //						}    
        //					}
        //				}
        ////#pragma omp barrier
        //				endtime = MPI_Wtime();
        //			//for(int num = 1; num <= (pWBR->maxIRow()-pWBR->minIRow()+1)*(pWBR->maxICol()-pWBR->minICol()+1); num++) 
        //			//{
        //			//	int iRow=num/(pWBR->maxICol()-pWBR->minICol()+1);
        //			//	int iCol=num%(pWBR->maxICol()-pWBR->minICol()+1);
        //			//	if (iCol==0)
        //			//	{
        //			//		iRow=iRow+pWBR->minIRow()-1;
        //			//		iCol=pWBR->maxICol();
        //			//	} 
        //			//	else
        //			//	{
        //			//		iRow=iRow+pWBR->minIRow();
        //			//		iCol=iCol+pWBR->minICol()-1;
        //			//	}
        //
        //			//	CellCoord coord(iRow, iCol);
        //			//	//if(!Operator(coord)) 
        //			//	//{
        //			//	//	cout<<"Operator is not sucessess!"<<endl;
        //			//	//	flag = false;
        //			//	//	break;
        //
        //			//	//}    
        //			//	Operator(coord);
        //			//}
        //	
        //			starttime2 = MPI_Wtime();
        //#pragma omp for 
        //			for(int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) 
        //			{
        //				for(int iCol = pWBR->minICol();	iCol <= pWBR->maxICol(); iCol++) 
        //				{
        //					CellCoord coord(iRow, iCol);
        //					if(!updatematrix(coord)) 
        //					{
        //						cout<<"updatematrix is not sucessess!"<<endl;
        //						flag = false;
        //						break;
        //
        //					}    
        //				}
        //			}
        //			}
        //
        //			endtime2 = MPI_Wtime();
        //			tcount+=endtime-starttime;
        //			tcount2+=endtime2-starttime2;
        //			COMNI.rowComm();
        //			MPI_Allreduce(&Termination, &noterm, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        //			
        //			itera++;
        //		}while( !noterm );
        //		MPI_Barrier(MPI_COMM_WORLD);
        //		cout<<"calculate matrix time:"<<tcount<<endl;
        //		cout<<"update matrix time:"<<tcount2<<endl;
        //		cout<<"iterative numerber is "<<itera<<endl;
    }
    else if (Application::_programType == CUDA_Type) {
    }
    else {
    }

    return flag;
}

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Run() {
    if (Work(_pWorkBR)) {
        return true;
    }
    return false;
}

#endif
