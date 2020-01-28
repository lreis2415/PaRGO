#ifndef ComputLayer_H
#define ComputLayer_H

/***************************************************************************
* computLayer.h
*
* Project: GPRO, v 2.0
* Purpose: Header file for class GPRO::ComputLayer
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2017. Ai Beibei
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include "rasterLayer.h"
#include "utility.h"
#include <iostream>

using namespace std;
#define Eps 0.0000001

inline bool ifDecomposeBySpace(const string& arg) {
    if (StringMatch(arg, "space")) {
        return true;
    }
    if (StringMatch(arg, "compute")) {
        return false;
    }
    return true;
}

namespace GPRO {
    template<class elemType>
    class ComputLayer : public RasterLayer<elemType> {
    public:
        ComputLayer();
		ComputLayer( const string RasterLayerName = "Untitled" );
		ComputLayer( RasterLayer<elemType> * dataLayers,
			         const int tmpGrain,
			         const string RasterLayerName = "Untitled" );
        ComputLayer( vector<RasterLayer<elemType> *> dataLayers,
                     const int tmpGrain,
                     const string RasterLayerName = "Untitled" );
        ~ComputLayer();

        int getGrain() { return _comptGrain; }
        void setComputGrain(int comptGrain) { _comptGrain = comptGrain;}
        void cleanDataLayers();
        vector<RasterLayer<elemType> *> *dataLayers();
        const vector<RasterLayer<elemType> *> *dataLayers() const;
        void addRasterLayer(RasterLayer<elemType> &dataLayer);

        bool newMetaData( int comptGrain );    //参数冗余，已在数据成员中添加
        bool getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR );

        bool readComptFile( const char *loadFile, const char* nbrFile);
        bool writeComptFile( const char *outputfile );
    public:
        vector<RasterLayer<elemType> *> _pDataLayers;
    protected:
        double _comptGrain;
    };
};

template<class elemType>
inline GPRO::ComputLayer<elemType>::
ComputLayer()
    :RasterLayer<elemType>() {
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
ComputLayer( const string RasterLayerName )
    :RasterLayer<elemType>( RasterLayerName ) {
    //_comptGrain = 10;
    //_pDataLayers.push_back(this);
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
	ComputLayer( RasterLayer<elemType> * dataLayer, const int tmpGrain, const string RasterLayerName )
	: RasterLayer<elemType>( RasterLayerName ), _comptGrain( tmpGrain ) {
		_pDataLayers.push_back(dataLayer);
    //GPRO::RasterLayer<elemType>::_pMetaData=dataLayer->metaData(); // wyj 12-19 太乱，临时加一笔看能不能跑。发现不能，存了人家的引用导致多重析构了
		//newMetaData(_comptGrain);	//没有邻域信息，故而这里无法new出workBR ///wyj 所以要再手动调用newMetaData...
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
ComputLayer( vector<RasterLayer<elemType> *> dataLayers, const int tmpGrain, const string RasterLayerName )
    : RasterLayer<elemType>( RasterLayerName ),
      _pDataLayers( dataLayers ), _comptGrain( tmpGrain ) {
    //newMetaData(_comptGrain);	//没有邻域信息，故而这里无法new出workBR ///wyj 所以要再手动调用newMetaData...
}

template<class elemType>
inline GPRO::ComputLayer<elemType>::
~ComputLayer() {
    //这里会自动调用基类的析构函数，只需要再释放_pDataLayers成员即可
    cleanDataLayers();
}

template<class elemType>
void GPRO::ComputLayer<elemType>::
cleanDataLayers() {
    //dataLayers中存放的是图层的指针，这里值释放vector，并不真正释放这些指针所指的图层
    //这些图层后续会继续使用，直到计算结束，调用自身析构函数去释放
    vector<RasterLayer<elemType> *> vTemp;
    vTemp.swap( _pDataLayers );
}

template<class elemType>
inline vector<GPRO::RasterLayer<elemType> *> *GPRO::ComputLayer<elemType>::
dataLayers() {
    return &_pDataLayers;
}
template<class elemType>
inline void GPRO::ComputLayer<elemType>::
addRasterLayer(RasterLayer<elemType> &dataLayer) {
    _pDataLayers.push_back(dataLayer);
}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
newMetaData( int comptGrain ) {
    //目前仅有进程0会来调用此函数
    //_pMetaData = new MetaData();	//VS这样的方式通过，而gc++不允许直接使用来自基类部分的数据成员
    if ( _pDataLayers.empty()) {
        return false;
    }
    const MetaData &rhs = *( _pDataLayers[0]->_pMetaData );

    RasterLayer<elemType>::_pMetaData = new MetaData();    //new元数据
    MetaData *&pMetaData = RasterLayer<elemType>::_pMetaData;    //指针的引用，只是为了简化书写,且用毕无需释放

    pMetaData->cellSize = rhs.cellSize * comptGrain;
    pMetaData->row = rhs._localworkBR.nRows() / comptGrain;    //这里的元数据都需要根据粒度换算
    pMetaData->row += ( rhs._localworkBR.nRows() % comptGrain ) ? 1 : 0;
    pMetaData->column = rhs._localworkBR.nCols() / comptGrain;
    pMetaData->column += ( rhs._localworkBR.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->format = rhs.format;
    pMetaData->projection = rhs.projection;
    pMetaData->noData = rhs.noData;
    pMetaData->myrank = rhs.myrank;
    pMetaData->processor_number = 0;    //目前只是串行构建
    pMetaData->_domDcmpType = rhs._domDcmpType;    //计算域的划分方式未必与数据域相同;目前是串行的,故是non_dcmp
    SpaceDims sdim( pMetaData->row, pMetaData->column );
    pMetaData->_glbDims = sdim;
    if ( pMetaData->_domDcmpType == NON_DCMP ) {
        CoordBR _glbWorkBR;
        RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims );    //根据计算域的邻域范围去求计算空间
        //计算域这里也只处理“数据范围-邻域范围”的范围
        pMetaData->_localworkBR = _glbWorkBR;
        CellCoord nwCorner( 0, 0 );
        CellCoord seCorner( pMetaData->_glbDims.nRows() - 1, pMetaData->_glbDims.nCols() - 1 );
        CoordBR subMBR( nwCorner, seCorner );
        pMetaData->_MBR = subMBR;
        pMetaData->_localdims = pMetaData->_glbDims;
    } else {
        cerr << "not support computLayer parallel construct now." << endl;
        return false;
    }

    pMetaData->dataType = RasterLayer<elemType>::getGDALType();

    for ( int i = 0; i < 6; i++ ) {
        pMetaData->pTransform[i] = rhs.pTransform[i];
    }
    pMetaData->pTransform[0] += rhs._localworkBR.minICol() * rhs.cellSize;//计算域左上角坐标是工作空间范围开始的
    pMetaData->pTransform[3] -= rhs._localworkBR.minIRow() * rhs.cellSize;
    pMetaData->pTransform[1] *= comptGrain;//东西、南北方向一个像素对应的距离，需更新
    pMetaData->pTransform[5] *= comptGrain;

    RasterLayer<elemType>::newCellSpace( pMetaData->_localdims, 0 ); //allocate,计算域栅格值初始化为0

    return true;
}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR ) {
    int myRank, process_nums;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &process_nums );

    int *pDcmpIdx = new int[nSubSpcs * 4];    //存储块划分结果，并广播
    if ( myRank == 0 ) {    //目前负载的划分是由主进程独立完成
        vector<CoordBR> vDcmpBR;    //所有划分快,共nSubSpcs个;因为MPI不支持通信自定义类型，故实际无用
        //根据邻域，主进程先读计算域所需图层数据，求解出工作空间范围（数据图层是有邻域成员的，计算域图层暂时没有或与其一致），再根据粒度，创建计算域图层；（总范围与总工作空间范围左上对齐，基本一致）
        //调用求解函数，即for循环compuLayer的每个栅格，配合粒度值，求解出compuLayer的各栅格值
        //创建decomp对象，调用划分函数valRowDcmp()，根据compuLayer值对其进行范围划分，结果由vector<CoordBR>& 返回
        //根据compuLayer图层的划分结果，映射到数据的工作空间范围，返回给主函数的vDcmpIdx
        //先串行求解，主进程更新了metedata后，通信给各子进程

        vector<CoordBR> vComptDcmpBR;
        cout << "hold " << RasterLayer<elemType>::_pMetaData->_glbDims << endl;
        DeComposition<elemType> deComp( RasterLayer<elemType>::_pMetaData->_glbDims, *( RasterLayer<elemType>::_pNbrhood ));
        if ( dcmpType == ROWWISE_DCMP ) {
            deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs );    //按值划分，故需要图层为参数;划分结果会以引用传回给vComptDcmpBR
            //_pDataLayers[0]->_pMetaData->_domDcmpType = ROWWISE_DCMP;	//暂时;估计没用//0316删除，不知道是否影响
        } else {
            cerr << "computLayer L388: not support until now." << endl;
        }
        //将划分结果映射给数据空间的子范围
        CoordBR glbWorkBR;
        Neighborhood<elemType> *pDataNbrhood = _pDataLayers[0]->nbrhood();
        pDataNbrhood->calcWorkBR( glbWorkBR, _pDataLayers[0]->_pMetaData->_glbDims );    //数据图层的全局工作空间
        int subBegin = glbWorkBR.minIRow(), subEnd = glbWorkBR.minIRow() - 1;
        int outputRows=_pDataLayers[0]->_pMetaData->_glbDims.nRows();
        int loadFileRows=this->metaData()->_glbDims.nRows();
        _comptGrain=outputRows/double(loadFileRows);//wyj 2019-12-6 加了这一行...解决读负载图时，图源粒度不统一的问题（比如空间均衡划分的负载图是等比例的，但计算负载划分的负载图是1:10的）
        int i = 0;
        for ( ; i < nSubSpcs - 1; ++i ) {
            subBegin = vComptDcmpBR[i].minIRow() * _comptGrain + glbWorkBR.minIRow();
            subEnd = ( vComptDcmpBR[i + 1].minIRow()) * _comptGrain + glbWorkBR.minIRow() - 1;
            //cout<<i<<" "<<subBegin<<" , "<<subEnd<<endl;
            CellCoord nwCorner( subBegin, glbWorkBR.minICol());
            CellCoord seCorner( subEnd, glbWorkBR.maxICol());
            CoordBR subMBR( nwCorner, seCorner );
            vDcmpBR.push_back( subMBR );
            pDcmpIdx[4 * i] = subBegin;
            pDcmpIdx[4 * i + 1] = glbWorkBR.minICol();
            pDcmpIdx[4 * i + 2] = subEnd;
            pDcmpIdx[4 * i + 3] = glbWorkBR.maxICol();
        }
        //存储转换，返回
        CellCoord nwCorner( subEnd + 1, glbWorkBR.minICol());
        CellCoord seCorner( glbWorkBR.maxIRow(), glbWorkBR.maxICol());
        CoordBR subMBR( nwCorner, seCorner );
        vDcmpBR.push_back( subMBR );
        pDcmpIdx[4 * i] = subEnd + 1;
        pDcmpIdx[4 * i + 1] = glbWorkBR.minICol();
        pDcmpIdx[4 * i + 2] = glbWorkBR.maxIRow();
        pDcmpIdx[4 * i + 3] = glbWorkBR.maxICol();
    }

    MPI_Barrier( MPI_COMM_WORLD );    //目前设计中，主进程会去完成计算域图层的计算，各进程在这里等候，计算完成后获取各自的subWorkBR即可
    //目前是串行，故需要广播，若各进程都来划分一次or并行划分，则需修改
    MPI_Bcast( pDcmpIdx, process_nums * 4, MPI_INT, 0, MPI_COMM_WORLD );
    CellCoord nwCorner2( pDcmpIdx[4 * myRank], pDcmpIdx[4 * myRank + 1] );
    CellCoord seCorner2( pDcmpIdx[4 * myRank + 2], pDcmpIdx[4 * myRank + 3] );
    CoordBR tmpWorkBR( nwCorner2, seCorner2 );
    subWorkBR = tmpWorkBR;    //已测试，这样赋值没问题
    //cout<<"tmpWorkBR "<<tmpWorkBR<<" subWorkBR "<<subWorkBR<<" shoule be the same"<<endl;
    MPI_Barrier( MPI_COMM_WORLD );
    delete[]pDcmpIdx;

    return true;
}
template<class elemType>
bool GPRO::ComputLayer<elemType>::
readComptFile(const char *loadFile, const char* nbrFile ) {
    RasterLayer<elemType> loadLayer("loadLayer");
    loadLayer.readNeighborhood(nbrFile);
    loadLayer.readFile(loadFile);
    loadLayer.newCellSpace(loadLayer.metaData()->_localdims);
    addRasterLayer(loadLayer);
    return true;
}

template<class elemType>
bool GPRO::ComputLayer<elemType>::
writeComptFile( const char *outputfile ) {
    //目前仅支持串行写出
    GDALAllRegister();

    if ( !RasterLayer<elemType>::createFile( outputfile )) {
        cout << "create file is not correct!" << endl;
        MPI_Finalize();
    }

    GDALDataset *poDataset = NULL;
    poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
    if ( poDataset == NULL /*检查是否正常打开文件*/) {
        //do something
        cout << "data file is not open correct" << endl;
        exit( 1 );
    }

    GDALRasterBand *poBanddest = poDataset->GetRasterBand( 1 );
    if ( poBanddest == NULL ) {
        //do something
        cout << "poBanddest is NULL" << endl;
        exit( 1 );
    }
    poBanddest->SetNoDataValue( RasterLayer<elemType>::_pMetaData->noData );

    if ( RasterLayer<elemType>::_pMetaData->myrank == 0 ) {
        poBanddest->RasterIO( GF_Write,
                              0,
                              0,
                              RasterLayer<elemType>::_pMetaData->_glbDims.nCols(),
                              RasterLayer<elemType>::_pMetaData->_glbDims.nRows(),
                              RasterLayer<elemType>::_pCellSpace->_matrix,
                              RasterLayer<elemType>::_pMetaData->_glbDims.nCols(),
                              RasterLayer<elemType>::_pMetaData->_glbDims.nRows(),
                              RasterLayer<elemType>::_pMetaData->dataType,
                              0,
                              0 );
    }

    if ( poDataset != NULL ) {
        GDALClose((GDALDatasetH) poDataset );
        poDataset = NULL;
    }

    return true;
}

#endif