#ifndef Transformation_H
#define Transformation_H

/***************************************************************************
* transformation.h
*
* Project: GPRO, v 2.0
* Purpose: Header file for class GPRO::Transformation
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/
#include <vector>
#include "basicTypes.h"
#include "cellSpace.h"
#include "application.h"
#include "computLayer.h"
#include <iostream>

using namespace std;

namespace GPRO
{
	template <class elemType>
	class Transformation
	{
	public:
		Transformation();
		Transformation( ComputLayer<elemType>* pLayer );	//继承类自定义实现时调用此类型
		Transformation( elemType load1, elemType load2, ComputLayer<elemType>* pLayer );
		//基类根据负载实现了一个默认版本，用于数据分布不均匀而计算分布均匀的情况

		virtual ~Transformation(){}	//可能是作为基类，因此写成虚析构;防止基类指针指向子类对象，释放不当

		virtual bool isTermination(bool isIter=false) {return isIter;}	//继承类直接改吧termination成员值
		bool Configure(ComputLayer<elemType>* pLayer, bool isCommunication);
		bool paramInit();
		virtual bool Operator(const CellCoord &coord);
		bool run();

	private:
		elemType minLoad;
		elemType maxLoad;
		int _myRank;

		double _noData;	//数据图层的空值
		int _computGrain;
		CoordBR _dataMBR;	//数据图层的行列号
		CoordBR _dataWorkBR;	//本进程处理的块对应的数据图层范围
		//int maxDataRow;
		//int maxDataCol;
		//int glbBeginRow;
		//int glbBeginCol;

	protected:
		//vector<RasterLayer<elemType>* > CommVec;	//待定，若为并行求解，则会需要
		ComputLayer<elemType>* _pComptLayer;
		//vector<RasterLayer<elemType>* > _pDataLayersV;	//引入这个是为了实现一些默认强度版本
		CoordBR* _pWorkBR;
		int Termination;
	};
};

template <class elemType>
inline GPRO::Transformation<elemType>::
Transformation()
	:minLoad(0),
	maxLoad(0),
	_pComptLayer(NULL),
	_pWorkBR(NULL),
	Termination(1)
{
}


template <class elemType>
inline GPRO::Transformation<elemType>::
Transformation( ComputLayer<elemType>* pLayer )	//继承类自定义实现时调用此类型
	:minLoad(0),
	maxLoad(0),
	_pComptLayer(pLayer),
	_pWorkBR(NULL),
	Termination(1)
{
	Configure(pLayer,false);
	_myRank = pLayer->_pMetaData->myrank;
}


template <class elemType>
inline GPRO::Transformation<elemType>::
Transformation( elemType load1, elemType load2, ComputLayer<elemType>* pLayer )	//基类根据负载指定默认版本
	:minLoad(load1),
	maxLoad(load2),
	_pComptLayer(pLayer),
	_pWorkBR(NULL),
	Termination(1)
{
	Configure(pLayer,false);
	_myRank = pLayer->_pMetaData->myrank;
}


template <class elemType>
bool GPRO::Transformation<elemType>::
Configure(ComputLayer<elemType>* pLayer, bool isCommunication)
{
	//目前很简略，以后可能会设计通信及更多数据成员
	if(_pWorkBR == NULL)
	{
		_pWorkBR = &pLayer->_pMetaData->_localworkBR;
	}
	if(isCommunication)
	{
		//do sth.
	}

	return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
paramInit()
{
	//必要的private数据成员一次性初始化，operator函数中会用到
	if( _pComptLayer->_pDataLayers.empty() ){
		cerr<<"Datalayers used for calculating compute layer should not be null."<<endl;
		return false;
	}
	_noData = _pComptLayer->_pDataLayers[0]->_pMetaData->noData;
	_computGrain = (int)(_pComptLayer->_pMetaData->cellSize / _pComptLayer->_pDataLayers[0]->_pMetaData->cellSize);
	_dataMBR = _pComptLayer->_pDataLayers[0]->_pMetaData->_MBR;
	_dataWorkBR = _pComptLayer->_pDataLayers[0]->_pMetaData->_localworkBR;

	return true;
}

//基类根据负载实现了一个默认版本，用于数据分布不均匀而计算分布均匀的情况;且计算分布是可指定的，不是复杂计算
template <class elemType>
bool GPRO::Transformation<elemType>::
Operator(const CellCoord &coord)
{
	//根据给定负载分布，进行计算
	int cRow = coord.iRow();
	int cCol = coord.iCol();
	CellSpace<elemType> &computL = *(_pComptLayer->cellSpace());

	if( cRow == _pWorkBR->minIRow() && cCol == _pWorkBR->minICol() ){
		cout<<"Transformation operator() function called."<<endl;
		//调用一个数据成员初始化函数，对行列范围等一次性初始化
		if( minLoad == maxLoad ){
			cerr<<"The load is balanced. No need to use compute layer."<<endl;
			return false;
		}
		if( minLoad<0 || maxLoad<0 ){
			cerr<<"The load specified cannot be negative."<<endl;
			return false;
		}

		if( !paramInit() ){
			return false;
		}
	}
	computL[cRow][cCol] = 0.0;	//初始化
	for( typename vector<RasterLayer<elemType>* >::iterator iter = _pComptLayer->_pDataLayers.begin(); iter!=_pComptLayer->_pDataLayers.end(); ++iter ){
		//对每个图层遍历计算，累积给计算域图层值
		CellSpace<elemType> &dataL = *((*iter)->cellSpace());	//模板类迭代指针这样用是否正确
		for( int dRow = cRow*_computGrain+_dataWorkBR.minIRow(); dRow<(cRow+1)*_computGrain+_dataWorkBR.minIRow(); ++dRow ){
			for( int dCol = cCol*_computGrain+_dataWorkBR.minICol(); dCol<(cCol+1)*_computGrain+_dataWorkBR.minICol(); ++dCol ){
				if( dRow > _dataMBR.maxIRow()||dRow>=dataL.nRows() || dCol > _dataMBR.maxICol() ||dCol>=dataL.nCols()){
					continue;
				}else{
					if( fabs(dataL[dRow][dCol] - _noData)>Eps && fabs(dataL[dRow][dCol] + 9999)>Eps){	//9999是针对我们的测试数据而多写的，其实没必要，属于数据问题，不属于程序问题
						computL[cRow][cCol] += maxLoad;//wyj 2019-11-12 noData不应该+=minLoad吗
					}else{
						computL[cRow][cCol] += minLoad;
					}
				}
			}
		}
	}

	return true;
}

template <class elemType>
bool GPRO::Transformation<elemType>::
run()
{
	if( Application::_programType != MPI_Type && Application::_programType != MPI_OpenMP_Type ){
		cerr<<"not supported yet."<<endl;
	}
	bool flag = true;
	//MPI_Barrier(MPI_COMM_WORLD);
	if( _myRank == 0 ){	//目前仅支持串行求解
		int iterNum = 0;	//迭代次数
		int termSum = 1;
		do
		{
			Termination = 1;
			for(int iRow = _pWorkBR->minIRow(); iRow <= _pWorkBR->maxIRow(); iRow++) 
			{
				for(int iCol = _pWorkBR->minICol();	iCol <= _pWorkBR->maxICol(); iCol++) 
				{
					CellCoord coord(iRow, iCol);
					if(!Operator(coord)) 
					{
						cout<<"Operator is not successes!"<<endl;
						flag = false;
						break;
					}
				}
			}
			//MPI_Allreduce(&Termination, &termSum, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
			termSum = Termination;	//暂提供给目前的串行版本
		} while (!termSum);
	}
	//MPI_Barrier(MPI_COMM_WORLD);

	return flag;
}


#endif
