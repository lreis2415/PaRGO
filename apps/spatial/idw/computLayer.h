#ifndef ComputLayer_H
#define ComputLayer_H

#include "rasterLayer.h"
#include <iostream>
using namespace std;
#define Eps 0.0000001

namespace GPRO{
	template <class elemType>
	class ComputLayer : public RasterLayer<elemType>
	{
	public:
		ComputLayer();
		ComputLayer( const string RasterLayerName = "Untitled");
		ComputLayer( vector<RasterLayer<elemType>* > dataLayers, const string RasterLayerName = "Untitled");
		~ComputLayer();

		void cleanDataLayers();
		vector<RasterLayer<elemType>* >* dataLayers();
		const vector<RasterLayer<elemType>* >* dataLayers() const;	//这里返回地址不可以是const,why
		//第一个const,修饰返回值；第二个const,类成员函数为const，即此函数不会修改数据成员

		bool newMetaData( const MetaData& rhs, int compuSize );
		bool transformation();
		bool getCompuLoad( int* pDcmpIdx, DomDcmpType dcmpType, int computSize, const int nSubSpcs );

		bool writeComptFile(const char* outputfile);

	protected:
		vector<RasterLayer<elemType>* > _pDataLayers;
	};
};

template <class elemType>
inline GPRO::ComputLayer<elemType>::
	ComputLayer()
	:RasterLayer<elemType>() 
{
}

template <class elemType>
inline GPRO::ComputLayer<elemType>::
	ComputLayer(const string RasterLayerName)
	:RasterLayer<elemType>( RasterLayerName )
{
}

template <class elemType>
inline GPRO::ComputLayer<elemType>::
	ComputLayer( vector<RasterLayer<elemType>* > dataLayers, const string RasterLayerName )
	: RasterLayer<elemType>( RasterLayerName ),
	_pDataLayers(dataLayers)
{
}

template <class elemType>
inline GPRO::ComputLayer<elemType>::
	~ComputLayer() 
{
	//这里会自动调用基类的析构函数，只需要再释放_pDataLayers成员即可
	cleanDataLayers();
}

template <class elemType>
void GPRO::ComputLayer<elemType>::
	cleanDataLayers() 
{
	//dataLayers中存放的是图层的指针，这里值释放vector，并不真正释放这些指针所指的图层
	//这些图层后续会继续使用，直到计算结束，调用自身析构函数去释放
	vector<RasterLayer<elemType>* > vTemp;
	vTemp.swap(_pDataLayers);
}

template <class elemType>
inline vector<GPRO::RasterLayer<elemType>* >* GPRO::ComputLayer<elemType>::
	dataLayers() 
{
	return &_pDataLayers;
}

//template <class elemType>
//inline vector<GPRO::RasterLayer<elemType>* >* GPRO::ComputLayer<elemType>::
//	dataLayers() const
//{
//	return &_pDataLayers;
//}

template <class elemType>
bool GPRO::ComputLayer<elemType>::
newMetaData( const MetaData& rhs, int compuSize ){
		//_pMetaData = new MetaData();	//VS这样的方式通过，而gc++不允许直接使用来自基类部分的数据成员
		RasterLayer<elemType>::_pMetaData = new MetaData();
		MetaData* &pMetaData = RasterLayer<elemType>::_pMetaData;	//指针的引用，只是为了简化书写,且用毕无需释放

		pMetaData->cellSize = rhs.cellSize * compuSize;
		pMetaData->row = rhs._localworkBR.nRows() / compuSize;	//这里的元数据都需要根据粒度换算
		pMetaData->row += (rhs._localworkBR.nRows() % compuSize) ? 1 : 0;
		pMetaData->column = rhs._localworkBR.nCols() / compuSize;
		pMetaData->column += (rhs._localworkBR.nCols() % compuSize) ? 1 : 0;
		pMetaData->format = rhs.format;
		pMetaData->projection = rhs.projection;
		pMetaData->noData = rhs.noData;
		pMetaData->myrank = rhs.myrank;
		//pMetaData->processor_number = rhs.processor_number;
		pMetaData->processor_number = 0;	//目前只是串行构建
		pMetaData->_domDcmpType = rhs._domDcmpType;	//计算域的划分方式未必与数据域相同;目前是串行的
		SpaceDims sdim(pMetaData->row, pMetaData->column);
		pMetaData->_glbDims = sdim;
		if( pMetaData->_domDcmpType == NON_DCMP ){
			CoordBR _glbWorkBR;
			//Neighborhood<elemType>* &pNbrhood = RasterLayer<elemType>::_pNbrhood;
			RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims );	//根据计算域的邻域范围去求计算空间
			//计算域这里也只处理“数据范围-邻域范围”的范围
			pMetaData->_localworkBR = _glbWorkBR;
			//cout<<"comptLayer L113 "<<*(RasterLayer<elemType>::_pNbrhood)<<endl;
			//int glbBegin = _glbWorkBR.nwCorner().iRow();
			//int glbEnd = _glbWorkBR.seCorner().iRow();
			//CellCoord nwCorner(glbBegin + pNbrhood->minIRow(),
			//	0);
			//CellCoord seCorner(glbEnd + pNbrhood->maxIRow(),
			//	pMetaData->_glbDims.nCols() - 1);
			CellCoord nwCorner(0, 0);
			CellCoord seCorner(pMetaData->_glbDims.nRows()-1, pMetaData->_glbDims.nCols()-1);
			CoordBR subMBR(nwCorner, seCorner);
			pMetaData->_MBR = subMBR;
			pMetaData->_localdims = pMetaData->_glbDims;
			//cout<<"comptLayer L127"<<" dcmpType "<<pMetaData->_domDcmpType<<" glbDims "<<pMetaData->_glbDims<<" localDims "<<pMetaData->_localdims<<" MBR "<<pMetaData->_MBR<<endl;
			//cout<<"comptLayer L128"<<" _localworkBR "<<pMetaData->_localworkBR<<endl;
		}else{
			cerr<<"not support computLayer decomposition now."<<endl;
			return false;
		}

		pMetaData->dataType = RasterLayer<elemType>::getType();

		for(int i = 0; i < 6; i++)
		{
			pMetaData->pTransform[i] = rhs.pTransform[i];
		}
		pMetaData->pTransform[0] += rhs._localworkBR.minICol()*rhs.cellSize;//计算域左上角坐标是工作空间范围开始的
		pMetaData->pTransform[3] -= rhs._localworkBR.minIRow()*rhs.cellSize;
		pMetaData->pTransform[1] *= compuSize;//东西、南北方向一个像素对应的距离，需更新
		pMetaData->pTransform[5] *= compuSize;

		//newCellSpace(pMetaData->_localdims,pMetaData->noData); //allocate
		RasterLayer<elemType>::newCellSpace(pMetaData->_localdims,0); //allocate,计算域栅格值初始化为0

		return true;
}

//用户接口函数，可考虑单写一个transformation类来继承，降低用户的已知信息量，尽量透明
template <class elemType>
bool GPRO::ComputLayer<elemType>::
transformation(){
	//这里写计算强度函数，求解计算域图层栅格值
	//cout<<RasterLayer<elemType>::_pMetaData->_localworkBR<<endl;
	CoordBR &workBR = RasterLayer<elemType>::_pMetaData->_localworkBR;
	//cout<<"L159 done"<<workBR.minIRow()<<" "<<workBR.maxIRow()<<" "<<workBR.minICol()<<" "<<workBR.maxICol()<<endl;
	CellSpace<elemType> &computL = *(RasterLayer<elemType>::_pCellSpace);
	CellSpace<elemType> &dataL = *(_pDataLayers[0]->cellSpace());
	double dataNoData = _pDataLayers[0]->_pMetaData->noData;
	int computSize = RasterLayer<elemType>::_pMetaData->cellSize / _pDataLayers[0]->_pMetaData->cellSize;
	//计算域图层不存在空值，只有0值
	//cout<<dataNoData<<endl;
	//dataNoData = -9999;
	int maxDRow = _pDataLayers[0]->_pMetaData->row;
	int maxDCol = _pDataLayers[0]->_pMetaData->column;
	int glbBeginRow = _pDataLayers[0]->_pMetaData->_localworkBR.minIRow();
	int glbBeginCol = _pDataLayers[0]->_pMetaData->_localworkBR.minICol();
	//cout<<"L171 done"<<maxDRow<<" "<<maxDCol<<" "<<glbBeginRow<<" "<<glbBeginCol<<endl;
	for(int cRow = workBR.minIRow(); cRow <= workBR.maxIRow(); cRow++) 
	{
		for(int cCol = workBR.minICol(); cCol <= workBR.maxICol(); cCol++) 
		{
			//cout<<computL[cRow][cCol]<<endl;
			//以本计算域cell内有效数据量为负载计算
			for( int dRow = cRow*computSize+glbBeginRow; dRow < (cRow+1)*computSize+glbBeginRow; ++dRow ){
				for( int dCol = cCol*computSize+glbBeginCol; dCol < (cCol+1)*computSize+glbBeginCol; ++dCol ){
					if( dRow > maxDRow-1 || dCol > maxDCol-1 ){
						continue;
					}else{
						if( fabs(dataL[dRow][dCol] - dataNoData)>Eps && fabs(dataL[dRow][dCol] + 9999)>Eps){	//9999是针对我们的测试数据而多写的，其实没必要，属于数据问题，不属于程序问题
							computL[cRow][cCol] += 10;
						}else{
							computL[cRow][cCol] += 1;
						}
					}
				}
			}
		}
	}
	return true;
}

template <class elemType>
bool GPRO::ComputLayer<elemType>::
getCompuLoad( int* pDcmpIdx, DomDcmpType dcmpType, const int computSize, const int nSubSpcs ){
		vector<CoordBR> vDcmpBR;
		//根据邻域，主进程先读计算域所需图层数据，求解出工作空间范围（数据图层是有邻域成员的，计算域图层暂时没有或与其一致），再根据粒度，创建计算域图层；（总范围与总工作空间范围左上对齐，基本一致）
		//调用求解函数，即for循环compuLayer的每个栅格，配合粒度值，求解出compuLayer的各栅格值
		//创建decomp对象，调用划分函数valRowDcmp()，根据compuLayer值对其进行范围划分，结果由vector<CoordBR>& 返回
		//根据compuLayer图层的划分结果，映射到数据的工作空间范围，返回给主函数的vDcmpIdx
		//先串行求解，主进程更新了metedata后，通信给各子进程
		if( _pDataLayers.empty() ){
			return false;
		}

		newMetaData( *(_pDataLayers[0]->_pMetaData), computSize );	//初始化了基类rasterLayer部分数据成员
		//cout<<"L207 done"<<RasterLayer<elemType>::_pMetaData->cellSize<<" "<<RasterLayer<elemType>::_pMetaData->row<<" "<<RasterLayer<elemType>::_pMetaData->column<<endl;
		//comptLayer覆盖的范围只有metaLayer的glbWorkBR范围
		transformation();	//求解computLayer._pCellSpace数据成员
		vector<CoordBR> vComptDcmpBR;	//待改名为vComptDcmpBR
		DeComposition<elemType> deComp(RasterLayer<elemType>::_pMetaData->_glbDims, *(RasterLayer<elemType>::_pNbrhood));
		if( dcmpType == ROWWISE_DCMP ){
			deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs);	//按值划分，故需要图层为参数;划分结果会以引用传回给vComptDcmpBR
			_pDataLayers[0]->_pMetaData->_domDcmpType = ROWWISE_DCMP;
		}else{
			cerr<<"not support until now."<<endl;
		}
		//for( vector<CoordBR>::iterator iter = vComptDcmpBR.begin(); iter!=vComptDcmpBR.end(); ++iter ){
		//	cout<<*iter<<endl;
		//}
		//将划分结果映射给数据空间的子范围
		CoordBR _glbWorkBR;
		Neighborhood<elemType> *pDataNbrhood = _pDataLayers[0]->nbrhood();
		pDataNbrhood->calcWorkBR( _glbWorkBR, _pDataLayers[0]->_pMetaData->_glbDims );	//数据图层的全局工作空间
		int subBegin = _glbWorkBR.minIRow(), subEnd = _glbWorkBR.minIRow()-1;
		int i = 0;
		for( ;i<nSubSpcs-1; ++i ){
			subBegin = vComptDcmpBR[i].minIRow()*computSize + _glbWorkBR.minIRow();
			subEnd = (vComptDcmpBR[i+1].minIRow())*computSize + _glbWorkBR.minIRow()-1;
			//cout<<i<<" "<<subBegin<<" , "<<subEnd<<endl;
			CellCoord nwCorner(subBegin, _glbWorkBR.minICol());
			CellCoord seCorner(subEnd, _glbWorkBR.maxICol());
			CoordBR subMBR(nwCorner, seCorner);
			vDcmpBR.push_back(subMBR);
			pDcmpIdx[4*i] = subBegin;
			pDcmpIdx[4*i+1] = _glbWorkBR.minICol();
			pDcmpIdx[4*i+2] = subEnd;
			pDcmpIdx[4*i+3] = _glbWorkBR.maxICol();
		}
		CellCoord nwCorner(subEnd+1, _glbWorkBR.minICol());
		CellCoord seCorner(_glbWorkBR.maxIRow(), _glbWorkBR.maxICol());
		CoordBR subMBR(nwCorner, seCorner);
		vDcmpBR.push_back(subMBR);
		pDcmpIdx[4*i] = subEnd+1;
		pDcmpIdx[4*i+1] = _glbWorkBR.minICol();
		pDcmpIdx[4*i+2] = _glbWorkBR.maxIRow();
		pDcmpIdx[4*i+3] = _glbWorkBR.maxICol();
		//for( vector<CoordBR>::iterator iter = vDcmpBR.begin(); iter != vDcmpBR.end(); ++iter ){
		//	cout<<*iter<<endl;	//it's ok here
		//}
		//因为目前MPI不支持通信自定义类型，故vDcmpIdx暂时实际上并没用

		////测试如何访问各元数据及栅格值
		//CoordBR* pWorkBR = &(_pMetaData->_localworkBR);
		//cout<<pWorkBR->minIRow()<<" "<<pWorkBR->maxIRow()<<" "<<pWorkBR->minICol()<<" "<<pWorkBR->maxICol()<<endl;
		//for(int iRow = pWorkBR->minIRow(); iRow <= pWorkBR->minIRow()+1; iRow++) 
		//{
		//	for(int iCol = pWorkBR->minICol();	iCol <= pWorkBR->minICol()+5; iCol++) 
		//	{
		//		CellSpace<double> &demL = *(_pCellSpace);
		//		cout<<demL[iRow][iCol]<<endl;
		//	}
		//}

		cout<<"getCompuLoad"<<endl;

		return true;
}

template <class elemType>
bool GPRO::ComputLayer<elemType>::
writeComptFile(const char* outputfile)
{
	//目前仅支持串行写出
	GDALAllRegister();

	if(!RasterLayer<elemType>::createFile(outputfile))
	{
		cout<<"create file is not correct!"<<endl;
		MPI_Finalize();
	}
	GDALDataset* poDataset = NULL;
	poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
	if( poDataset == NULL /*检查是否正常打开文件*/)
	{
		//do something
		cout<<"data file is not open correct"<<endl;
		exit(1);
	}
	GDALRasterBand*	poBanddest = poDataset->GetRasterBand(1);
	if (poBanddest == NULL)
	{
		//do something
		cout<<"poBanddest is NULL"<<endl;
		exit(1);
	}
	poBanddest->SetNoDataValue(RasterLayer<elemType>::_pMetaData->noData);

	if (RasterLayer<elemType>::_pMetaData->myrank == 0)
		poBanddest->RasterIO(GF_Write, 0, 0, RasterLayer<elemType>::_pMetaData->_glbDims.nCols(),RasterLayer<elemType>::_pMetaData->_glbDims.nRows(), RasterLayer<elemType>::_pCellSpace->_matrix, RasterLayer<elemType>::_pMetaData->_glbDims.nCols(), RasterLayer<elemType>::_pMetaData->_glbDims.nRows(), RasterLayer<elemType>::_pMetaData->dataType, 0, 0);

	if (poDataset != NULL)
	{
		GDALClose((GDALDatasetH)poDataset);
		poDataset = NULL;
	}

	return true;
}


#endif