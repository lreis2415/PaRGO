/**
 * \file basicTypes
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::RasterOperator
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
#ifndef RASTEROPERATOR_H
#define RASTEROPERATOR_H

#include "basicTypes.h"
#include "basicCell.h"
#include "metaData.h"
#include "rasterLayer.h"
#include "communication.h"
#include "application.h"
#include <string>
#include <iostream>
#include <typeinfo>
#include <vector>
#include "mpi.h"
#include <omp.h>
#include <gdal_priv.h>

using namespace std;

namespace GPRO
{
    /**
     * \ingroup gpro
     * \class RasterOperator
     * \brief A basic super class that each Operator of a specific algorithm should extends
     */
    template <class elemType>
    class RasterOperator
    {
    public:
        RasterOperator()
            : _domDcmpType(NON_DCMP),
              _pCellSpace(NULL),
              _pNbrhood(NULL),
              _pWorkBR(NULL),
              _pComptLayer(NULL),
              commFlag(false),
              Termination(true) {
        }

        virtual ~RasterOperator() {
        }

        //virtual void processing(){}
        //virtual bool updatematrix(const CellCoord &coord){}

        /**
         * \brief Basic function in which serial-style algorithm should be writen
         * \param[in] coord coordinate of the central cell(point)
         * \param[in] operFlag an unused controlling variable.
         */
        virtual bool Operator(const CellCoord& coord, bool operFlag) { return operFlag; }

        virtual bool isTermination() { return false; }

        bool Work(const CoordBR* const pWorkBR);
        bool Run();
        bool Configure(RasterLayer<elemType>* pLayer, bool isCommunication);
        void comptLayer(RasterLayer<elemType>& layerD);


    private:
        DomDcmpType _domDcmpType;

    public:
        RasterLayer<elemType>* _pComptLayer; ///暂时捕捉真实计算时间用
        vector<RasterLayer<elemType> *> CommVec;
        CellSpace<elemType>* _pCellSpace;
        Neighborhood<elemType>* _pNbrhood; //目前这两个成员变量并没有使用
        CoordBR* _pWorkBR;
        bool commFlag;
        int Termination;
    };
};

template <class elemType>
void GPRO::RasterOperator<elemType>::
comptLayer(RasterLayer<elemType>& layerD) {
    _pComptLayer = &layerD;
    _pComptLayer->cellSpace()->initVals(0); //wyj 2019-12-7 这样的话fcm idw Operator里就不用每个初始化了
    Configure(_pComptLayer, false);
}

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Configure(RasterLayer<elemType>* pLayer, bool isCommunication) {
    if (_pWorkBR == NULL) {
        _pWorkBR = &pLayer->_pMetaData->_localworkBR;
        _domDcmpType = pLayer->_pMetaData->_domDcmpType;
    }
    //cout<<_pWorkBR->minIRow()<<" "<<_pWorkBR->minICol()<<" "<<_pWorkBR->maxIRow()<<" "<<_pWorkBR->maxICol()<<endl;
    if (isCommunication) {
        CommVec.push_back(pLayer);
        //cout<<"_pCommVec->size is "<<endl;
        if (commFlag == false) {
            commFlag = true;
        }
    }

    return true;
}

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Work(const CoordBR* const pWBR) {
    bool flag = true; //标识本函数是否得以正确执行，会返回给run函数

    //cout<<"pWBR->minIRow() "<<pWBR->minIRow()<<"  pWBR->maxIRow()"<<pWBR->maxIRow()<<endl;
    //cout<<"_pCommVec->size() "<<CommVec.size()<<endl;
    Communication<elemType> COMNI(&CommVec);
    if (commFlag) {
        COMNI.setBuffer();
    }
    int noterm = 1; //是否继续迭代，考虑换别的变量
    int itera = 0; //迭代次数
    int myRank;
    int nRow=pWBR->maxIRow()-pWBR->minIRow();
    int delim=20;
    int lastRowInterval=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    double startTime=MPI_Wtime();
    double endTime;
    double iterStartTime;
    double rowStartTime=MPI_Wtime();
    if (Application::_programType == MPI_Type) {
        MPI_Barrier(MPI_COMM_WORLD);
        do {
            iterStartTime=MPI_Wtime();
            Termination = 1;
            for (int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) {
                
                for (int iCol = pWBR->minICol(); iCol <= pWBR->maxICol(); iCol++) {
                    CellCoord coord(iRow, iCol);
                    if (!Operator(coord, true)) {
                        cout << "Operator is not successes!" << endl;
                        flag = false;
                        break;
                    }
                }
                
                int rowInterval= (iRow-pWBR->minIRow())/(double(nRow)/delim);
                if(rowInterval != lastRowInterval) {
                    lastRowInterval=rowInterval;
                    cout<<"rank"<<myRank<<"\t[";
                    for (int i = 0; i < lastRowInterval; ++i) {
                        cout<<".";
                    }
                    for (int i = 0; i < delim-lastRowInterval; ++i) {
                        cout<<" ";
                    }
                    endTime=MPI_Wtime();
                    cout<<"]"<<endTime-rowStartTime<<"s ("<<iRow-nRow/delim<<"~"<<iRow<<")"<<endl;
                    rowStartTime=MPI_Wtime();
                }
            }
            endTime=MPI_Wtime();
            cout<<"rank"<<myRank<<" iter time "<<endTime-iterStartTime<<"s"<<endl;
            startTime=MPI_Wtime();
            if (commFlag) {
                //这里应该是统一调用函数，在函数内部，根据元数据获取的划分方式再确定调用哪一个通信
                if (_domDcmpType == ROWWISE_DCMP) {
                    COMNI.rowComm();
                }
                if (_domDcmpType == COLWISE_DCMP) {
                    COMNI.colComm();
                }
            }
            MPI_Allreduce(&Termination, &noterm, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
            itera++;
            
            endTime=MPI_Wtime();
            cout<<"rank"<<myRank<<" commu time "<<endTime-startTime<<"s"<<endl;
        }
        while (!noterm);
        MPI_Barrier(MPI_COMM_WORLD);
        //cout<<"iterative number is "<<itera<<endl;
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
        ;
    }
    else {
        ;
    }

    return flag;
}

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Run() {
    if (Work(_pWorkBR)) {
        return true;
    }
    else {
        return false;
    }
}

#endif
