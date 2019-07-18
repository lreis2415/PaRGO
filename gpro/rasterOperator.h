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

#include "utility.h"
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
        :_pCellSpace(nullptr),
         _pNbrhood(nullptr),
		 pWorkBR(nullptr),
		 comm(false),
		 Termination(true)
         {}

      virtual ~RasterOperator() DEFAULT;

      virtual bool Operator(const CellCoord &coord,bool operFlag) {return operFlag;}

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
	  int Termination;;
  };
};

template <class elemType>
bool GPRO::RasterOperator<elemType>::
Configure(RasterLayer<elemType>* pLayer, bool isCommunication)
{
	if(pWorkBR == nullptr)
	{
		pWorkBR = &pLayer->_pMetaData->_localworkBR;
	}
	if(isCommunication)
	{
		CommVec.push_back( pLayer );
		//cout<<"_pCommVec->size is "<<endl;
		if(!comm)
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
	COMNI.setBuffer();
	//cout<<"COMNI is OK!!!!!!"<<endl;
	bool flag = true;
	int noterm = 1;
	int itera = 0;
	if(Application::_programType == MPI_Type)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		double endtime1;
		double endtime2;
		double starttime;
	
		do
		{
			//starttime = MPI_Wtime();
			Termination = 1;
			for(int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) 
			{
				for(int iCol = pWBR->minICol();	iCol <= pWBR->maxICol(); iCol++) 
				{
					CellCoord coord(iRow, iCol);
					if(!Operator(coord,true)) 
					{
						cout<<"Operator is not successes!"<<endl;
						flag = false;
						break;
					}    
				}
			}
			
			COMNI.rowComm();

			MPI_Allreduce(&Termination, &noterm, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
			
			itera++;
			/*if(itera%1000 == 0)
			{
			  cout<<"itera is "<<itera<<endl;
			}*/

		}while( !noterm );
		MPI_Barrier(MPI_COMM_WORLD);
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        if( myrank==0 )
		cout<<"iterative number is "<<itera<<endl;
	}
	else if (Application::_programType == MPI_OpenMP_Type)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		
		do
		{
			Termination = 1;
			#pragma omp parallel for
			for(int iRow = pWBR->minIRow(); iRow <= pWBR->maxIRow(); iRow++) 
			{
				for(int iCol = pWBR->minICol();	iCol <= pWBR->maxICol(); iCol++) 
				{
					CellCoord coord(iRow, iCol);
					if(!Operator(coord,true)) 
					{
						cout<<"Operator is not successes!"<<endl;
						flag = false;
						break;
					}    
				}
			}

				
			COMNI.rowComm();

			MPI_Allreduce(&Termination, &noterm, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
			
			itera++;

		}while( !noterm );
		MPI_Barrier(MPI_COMM_WORLD);
		cout<<"iterative number is "<<itera<<endl;
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