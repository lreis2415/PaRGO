#ifndef DeComposition_H
#define DeComposition_H

/***************************************************************************
* DeComposition.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::DeComposition
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
#include "cellSpace.h"
#include "neighborhood.h"
#include "metaData.h"

#define Eps 0.0000001

namespace GPRO {
    template<class elemType>
    class ComputLayer;

    template<class elemType>
    class DeComposition {
    public:
        DeComposition()
            : _pNbrhood( 0 ) {}

        DeComposition( const SpaceDims &cellSpaceDims, const Neighborhood <elemType> &nbrhood );

        ~DeComposition() {}

        bool rowDcmp( MetaData &metaData, int nSubSpcs ) const;
        bool colDcmp( MetaData &metaData, int nSubSpcs ) const;
        /*bool blockDcmp(MetaData &metaData,
                       int nRowSubSpcs,
                       int nColSubSpcs) const;*/
		//这个应该改成_val1DDcmp这样的形式，提高通用性
        bool valRowDcmp( vector<CoordBR> &vDcmpIdx, ComputLayer<elemType> &computLayer, int nSubSpcs ) const;

    private:
        void _smpl1DDcmp( int &subBegin, int &subEnd,
                          int glbBegin, int glbEnd,
                          int nSubSpcs, int iSubSpc ) const;

    private:
        SpaceDims _glbDims;
        CoordBR _glbWorkBR;
        const Neighborhood <elemType> *_pNbrhood;
    };
};

/****************************************
*             Private Methods           *
*****************************************/

template<class elemType>
void GPRO::DeComposition<elemType>::
_smpl1DDcmp( int &subBegin, int &subEnd,
             int glbBegin, int glbEnd,
             int nSubSpcs, int iSubSpc ) const {
    int glbRange = glbEnd - glbBegin + 1;
    int remainder = glbRange % nSubSpcs;
    int blockSize;
    if ( remainder == 0 ) {
        blockSize = glbRange / nSubSpcs;
        subBegin = glbBegin + iSubSpc * blockSize;
        subEnd = subBegin + blockSize - 1;
    } else {
        blockSize = glbRange / nSubSpcs + 1;
        if ( iSubSpc < remainder ) {
            subBegin = glbBegin + iSubSpc * blockSize;
            subEnd = subBegin + blockSize - 1;
        } else {
            if ( iSubSpc == remainder ) {
                subBegin = glbBegin + iSubSpc * blockSize;
                subEnd = subBegin + blockSize - 2;
            } else {
                subBegin = glbBegin
                    + remainder * blockSize
                    + ( iSubSpc - remainder ) * ( blockSize - 1 );
                subEnd = subBegin + blockSize - 2;
            }
        }
    }
}

/****************************************
*             Public Methods            *
*****************************************/

template<class elemType>
GPRO::DeComposition<elemType>::
DeComposition( const SpaceDims &cellSpaceDims, const Neighborhood <elemType> &nbrhood ) {
    if ( !nbrhood.calcWorkBR( _glbWorkBR, cellSpaceDims )) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: invalid Neighborhood to construct a DeComposition" \
			 << endl;
        exit( 1 );
    }
    _glbDims = cellSpaceDims;
    _pNbrhood = &nbrhood;
}


template<class elemType>
bool GPRO::DeComposition<elemType>::
rowDcmp( MetaData &metaData, int nSubSpcs ) const {
    if ( nSubSpcs < 1 || nSubSpcs > _glbWorkBR.nRows()) {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: invalid number of SubSpaces (" << nSubSpcs << ")" \
             << " considering the global working bounding rectangle (" << _glbWorkBR << ")" \
             << endl;
        return false;
    }

    //DomDcmpType dcmpType = ROWWISE_DCMP;
    //metaData._domDcmpType = dcmpType;

    int glbBegin = _glbWorkBR.nwCorner().iRow();
    int glbEnd = _glbWorkBR.seCorner().iRow();
    int iSubSpc = metaData.myrank;
    int subBegin, subEnd;
    _smpl1DDcmp( subBegin, subEnd,
                 glbBegin, glbEnd,
                 nSubSpcs, iSubSpc );

    CellCoord nwCorner( subBegin + _pNbrhood->minIRow(), 0 );
    CellCoord seCorner( subEnd + _pNbrhood->maxIRow(), _glbDims.nCols() - 1 );
    CoordBR subMBR( nwCorner, seCorner );
    metaData._MBR = subMBR;

    SpaceDims dims( subMBR.nRows(), subMBR.nCols() );
    metaData._localdims = dims;
    CoordBR workBR;
    if ( !_pNbrhood->calcWorkBR( workBR, dims )) {
        return false;
    }
    metaData._localworkBR = workBR;

    return true;
}

template<class elemType>
bool GPRO::DeComposition<elemType>::
colDcmp( MetaData &metaData, int nSubSpcs ) const {
    if ( nSubSpcs < 1 || nSubSpcs > _glbWorkBR.nCols()) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: invalid number of SubSpaces (" << nSubSpcs << ")" \
			 << " considering the global working bounding rectangle (" << _glbWorkBR << ")" \
			 << endl;
        return false;
    }

	//DomDcmpType dcmpType = COLWISE_DCMP;	//同上，是否需要
    //metaData._domDcmpType = dcmpType;

    int glbBegin = _glbWorkBR.nwCorner().iCol();
    int glbEnd = _glbWorkBR.seCorner().iCol();
    int iSubSpc = metaData.myrank;
    int subBegin, subEnd;
    _smpl1DDcmp( subBegin, subEnd,
                 glbBegin, glbEnd,
                 nSubSpcs, iSubSpc );

	CellCoord nwCorner( 0, subBegin + _pNbrhood->minICol());
    CellCoord seCorner( _glbDims.nRows() - 1, subEnd + _pNbrhood->maxICol());
    CoordBR subMBR( nwCorner, seCorner );
    metaData._MBR = subMBR;

    SpaceDims dims( subMBR.nRows(), subMBR.nCols());
    metaData._localdims = dims;
    CoordBR workBR;
    if ( !_pNbrhood->calcWorkBR( workBR, dims )) {
        return false;
    }
    metaData._localworkBR = workBR;
	//cout<<workBR.minIRow()<<" "<<workBR.minICol()<<" "<<workBR.maxIRow()<<" "<<workBR.maxICol()<<endl;

    return true;
}

//template <class elemType>
//bool GPRO::DeComposition<elemType>::
//blockDcmp(SubInfoVect &vpSubSpcInfos,
//          int nRowSubSpcs,
//          int nColSubSpcs) const {
//  vpSubSpcInfos.clear();
//  int nSubSpcs = nRowSubSpcs * nColSubSpcs;
//  if(nSubSpcs < 1 ||
//     nRowSubSpcs > _glbWorkBR.nRows() ||
//     nColSubSpcs > _glbWorkBR.nCols()) {
//    cerr << __FILE__ << " " << __FUNCTION__ \
//         << " Error: invalid number of subspaces" \
//         << endl;
//    return false;
//  }
//
//  DomDcmpType dcmpType = BLOCK_DCMP;
//  int glbRowBegin = _glbWorkBR.nwCorner().iRow();
//  int glbRowEnd = _glbWorkBR.seCorner().iRow();
//  int glbColBegin = _glbWorkBR.nwCorner().iCol();
//  int glbColEnd = _glbWorkBR.seCorner().iCol();
//  
//  for(int iRowSubSpc = 0; iRowSubSpc < nRowSubSpcs; iRowSubSpc++) {
//    int subRowBegin, subRowEnd;
//    _smpl1DDcmp(subRowBegin, subRowEnd,
//                glbRowBegin, glbRowEnd,
//                nRowSubSpcs, iRowSubSpc);
//    for(int iColSubSpc = 0; iColSubSpc < nColSubSpcs; iColSubSpc++) {
//      int subColBegin, subColEnd;
//      _smpl1DDcmp(subColBegin, subColEnd,
//                  glbColBegin, glbColEnd,
//                  nColSubSpcs, iColSubSpc);
//
//      CellCoord nwCorner(subRowBegin + _pNbrhood->minIRow(), 
//                         subColBegin + _pNbrhood->minICol());
//      CellCoord seCorner(subRowEnd + _pNbrhood->maxIRow(), 
//                         subColEnd + _pNbrhood->maxICol());
//
//      CoordBR subMBR(nwCorner, seCorner);
//      SpaceDims dims(subMBR.nRows(), subMBR.nCols());
//      CoordBR workBR;
//      if(!_pNbrhood->calcWorkBR(workBR, dims)) {
//        return false; 
//      }
//
//      vector<IntVect> mNbrSpcIDs(8);
//      for(int iDir = 0; iDir < 8; iDir++) {
//        int iRowNbrSpc = iRowSubSpc + *(nbrLocations + 2*iDir);
//        int iColNbrSpc = iColSubSpc + *(nbrLocations + 2*iDir + 1);
//        if(iRowNbrSpc >= 0 && iRowNbrSpc < nRowSubSpcs &&
//           iColNbrSpc >= 0 && iColNbrSpc < nColSubSpcs) {
//          int nbrSpcID = iRowNbrSpc * nColSubSpcs + iColNbrSpc;
//          mNbrSpcIDs[iDir].push_back(nbrSpcID);
//        }
//      }
//
//      int iSubSpc = iRowSubSpc * nColSubSpcs + iColSubSpc;
//      int workload = workBR.size();
//      SubInfoVect::iterator iInsert = 
//        upper_bound(vpSubSpcInfos.begin(),
//                    vpSubSpcInfos.end(),
//                    InfoWLPair(0, workload),
//                    greater<InfoWLPair>());
//      vpSubSpcInfos.insert(iInsert,
//        InfoWLPair(new CellSpaceInfo(iSubSpc, dcmpType, _glbDims,
//                                    subMBR, workBR, mNbrSpcIDs),
//                   workload));
//    } // End of iColSubSpc loop
//  } // End of iRowSubSpc loop
//  return true;
//}

template<class elemType>
bool GPRO::DeComposition<elemType>::
valRowDcmp( vector<CoordBR> &vDcmpBR, ComputLayer<elemType> &layer, int nSubSpcs ) const {
    //目前仅主进程会访问这个函数，串行对全区划分；划分的对象是layer的全区范围
    //若设计并行策略,则划分边界等都有待更新
    if ( nSubSpcs < 1 || nSubSpcs > _glbWorkBR.nRows()) {
        cerr << __FILE__ << " " << __FUNCTION__ \
			 << " Error: invalid number of SubSpaces (" << nSubSpcs << ")" \
			 << " considering the global working bounding rectangle (" << _glbWorkBR << ")" \
			 << endl;
        return false;
    }

    //根据计算域图层值，按行均匀划分，范围索引返回到vDcmpBR;基于负载的划分是针对glbdims来分，而不是localWorkBR
    CellSpace<elemType> &comptL = *( layer.cellSpace());
    double *rowComptLoad = new double[layer._pMetaData->_glbDims.nRows()];
    double totalComptLoad = 0.0;
    //这里目前只针对串行
    for ( int cRow = 0; cRow < layer._pMetaData->_glbDims.nRows(); cRow++ ) {
        rowComptLoad[cRow] = 0.0;
        for ( int cCol = 0; cCol < layer._pMetaData->_glbDims.nCols(); cCol++ ) {
			//这里为什么对-9999需要单独写出来，nodata会无效？
            if (( comptL[cRow][cCol] + 9999 ) > Eps && ( comptL[cRow][cCol] - layer._pMetaData->noData ) > Eps ) {
                rowComptLoad[cRow] += comptL[cRow][cCol];
            }
        }
        totalComptLoad += rowComptLoad[cRow];
    }
    //cout << layer._pMetaData->_glbDims.nRows() << " rows with totalComptLoad " << totalComptLoad << endl;

    double averComptLoad = totalComptLoad / ( nSubSpcs );
    double accuComptLoad = 0;
    //int subBegin = pWorkBR->minIRow();
    int subBegin = 0;
    int subEnd = subBegin;    //至少分配一行；是计算域的一行
    for ( int i = 0; i < nSubSpcs - 1; ++i ) {    //共寻找nSubSpcs-1个划分位置,每个划分位置都在剩下的局部最优
        double subComptLoad = 0.0;
        double minComptDiff = 0.0;
        for ( int cRow = subBegin; cRow <= layer._pMetaData->_MBR.maxIRow(); cRow++ ) {
            int rowGap = cRow - 0;
            subComptLoad += rowComptLoad[rowGap];
            accuComptLoad += rowComptLoad[rowGap];
            if ( subComptLoad <= averComptLoad ) {
                if ( fabs( minComptDiff ) < Eps || ( averComptLoad - subComptLoad ) < minComptDiff ) {
                    minComptDiff = averComptLoad - subComptLoad;
                }
            } else {
                if ( fabs( minComptDiff ) < Eps || ( subComptLoad - averComptLoad ) <= minComptDiff ) {
                    subEnd = cRow;
                    CellCoord nwCorner( subBegin, 0 );
                    CellCoord seCorner( subEnd, layer._pMetaData->_MBR.maxICol());
                    CoordBR subMBR( nwCorner, seCorner );
                    vDcmpBR.push_back( subMBR );
                    averComptLoad = ( totalComptLoad - accuComptLoad ) / ( nSubSpcs - i - 1 );
                    //cout<<"averComptLoad "<<averComptLoad*10000/totalComptLoad<<endl;
                    //cout<<"subComptLoad "<<subComptLoad<<endl;
                    subBegin = cRow + 1;
                    //subComptLoad = 0;
                } else {
                    subEnd = cRow - 1;
                    accuComptLoad -= rowComptLoad[rowGap];
                    CellCoord nwCorner( subBegin, 0 );
                    CellCoord seCorner( subEnd, layer._pMetaData->_MBR.maxICol());
                    CoordBR subMBR( nwCorner, seCorner );
                    vDcmpBR.push_back( subMBR );
                    subComptLoad -= rowComptLoad[rowGap];
                    averComptLoad = ( totalComptLoad - accuComptLoad ) / ( nSubSpcs - i - 1 );
                    //cout<<"averComptLoad "<<averComptLoad*10000/totalComptLoad<<endl;
                    subBegin = cRow;
                    //subComptLoad = 0;
                }
                //cout<<"the "<<rowGap<<" subComptLoad "<<subComptLoad<<endl;
                break;
            }
        }
    }
    //最后一块，从subBegin到最后一行，都是剩下的
    CellCoord nwCorner( subBegin, 0 );
    CellCoord seCorner( layer._pMetaData->_MBR.maxIRow(), layer._pMetaData->_MBR.maxICol());
    CoordBR subMBR( nwCorner, seCorner );
    vDcmpBR.push_back( subMBR );

    delete[]rowComptLoad;
    //for( vector<CoordBR>::iterator iter = vDcmpBR.begin(); iter!=vDcmpBR.end(); ++iter ){
    //	cout<<*iter<<endl;
    //}
    return true;
}

#endif
