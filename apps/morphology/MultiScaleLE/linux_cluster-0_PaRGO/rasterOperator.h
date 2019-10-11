#ifndef RASTEROPERATOR_H
#define RASTEROPERATOR_H

/***************************************************************************
* rasterOperator.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::RasterOperator
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
#include "rasterLayer.h"
#include "communication.h"
#include "application.h"
#include <string>
#include <iostream>
#include <typeinfo>
#include <vector>
#include"mpi.h"
#include<omp.h>
#include<gdal_priv.h>


using namespace std;

namespace GPRO 
{
  
  template <class elemType>
  class RasterOperator 
  {
    public:
      RasterOperator()
        :_pCellSpace(NULL),
         _pNbrhood(NULL),
		 pWorkBR(NULL),
		 comm(false),
		 Termination(true)
         {}

      virtual ~RasterOperator() {}

      virtual bool Operator(const CellCoord &coord) {}

	  virtual bool isTermination() {return false;}

      bool Work(const CoordBR *const pWorkBR);
	  bool Run();
	  bool Configure(RasterLayer<elemType>* pLayer, bool isCommunication);

    public:

	  vector<RasterLayer<elemType>* > CommVec;
      CellSpace<elemType> *_pCellSpace;
      Neighborhood<elemType> *_pNbrhood;
	  CoordBR* pWorkBR;
	  bool comm;
	  bool Termination;;
  };
};

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Configure(RasterLayer<elemType>* pLayer, bool isCommunication)
{
	if(pWorkBR == NULL)
	{
		pWorkBR = &pLayer->_pMetaData->_localworkBR;
	}
	if(isCommunication)
	{
		CommVec.push_back( pLayer );
		//cout<<"_pCommVec->size is "<<endl;
		if(comm == false)
		{
			comm = true;
		}
	}

	return true;
}

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Work(const CoordBR *const pWBR)
{
	//cout<<"pWBR->minIRow() "<<pWBR->minIRow()<<"  pWBR->maxIRow()"<<pWBR->maxIRow()<<endl;
	
	
	//cout<<"_pCommVec->size() "<<CommVec.size()<<endl;
	//cout<<"Operator is OK!"<<endl;
	Communication<elemType> COMNI(&CommVec);
	//cout<<"COMNI is OK!!!!!!"<<endl;
	bool flag = true;
	bool noterm = true;
	int itera = 0;
	if(Application::_programType == MPI_Type)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		do
		{
			//Termination = true;
			for(int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) 
			{
				for(int iCol = pWBR->minICol();	iCol <= pWBR->maxICol(); iCol++) 
				{
					CellCoord coord(iRow, iCol);
					if(!Operator(coord)) 
					{
						cout<<"Operator is not sucessess!"<<endl;
						flag = false;
						break;
					}    
				}
			}

			//cout<<"iterator is done one!!!!!!"<<endl;

			COMNI.rowComm();
			//bool noterm;
			//MPI_Allreduce(&Termination, &noterm, 1, MPI_CHAR, MPI_LAND, MPI_COMM_WORLD);

			//itera++;
		}while( !noterm );
		MPI_Barrier(MPI_COMM_WORLD);
		//cout<<"iterative numerber is "<<itera<<endl;
	}
	else if (Application::_programType == MPI_OpenMP_Type)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		do
		{
			#pragma omp parallel for
			for(int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) 
			{
				for(int iCol = pWBR->minICol();	iCol <= pWBR->maxICol(); iCol++) 
				{
					CellCoord coord(iRow, iCol);
					if(!Operator(coord)) 
					{
						cout<<"Operator is not sucessess!"<<endl;
						flag = false;
						break;
					}    
				}
			}

			COMNI.rowComm();

		}while( isTermination() );
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else if(Application::_programType == CUDA_Type)
	{
		;
	}
	else
	{
		;
	}

	return flag;
}


template <class elemType>
bool GPRO::RasterOperator<elemType>::
Run()
{
	//cout<<"begin run!"<<endl;
	if(Work(pWorkBR))
	{
		//cout<<"Work is done!"<<endl;
		return true;
	}
	else
	{
		return false;
	}

}




#endif