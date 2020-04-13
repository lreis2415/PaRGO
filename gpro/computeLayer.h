/**
 * \file computLayer
 * \author Ai Beibei, Wang Yujing
 * \brief Header file for class GPRO::ComputeLayer
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

#ifndef ComputeLayer_H
#define ComputeLayer_H

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
    /**
     * \ingroup gpro
     * \class ComputeLayer
     * \brief the compute layer, stores the intensity values
     */
    template<class elemType>
    class ComputeLayer : public RasterLayer<elemType> {
    public:
        ComputeLayer();
		ComputeLayer( const string RasterLayerName = "Untitled" );
		ComputeLayer( RasterLayer<elemType> * dataLayers,
			         const int tmpGrain,
			         const string RasterLayerName = "Untitled" );
        ComputeLayer( vector<RasterLayer<elemType> *> dataLayers,
                     const int tmpGrain,
                     const string RasterLayerName = "Untitled" );
        
        ~ComputeLayer();

        /*
        * getters and setters
        */
        int getGrain() { return _comptGrain; }
        void setComputGrain(int comptGrain) { _comptGrain = comptGrain;}
        
        /**
         * \brief Init this compute layer. It should be called before solving the compute layer.
         *  
         * Methods with suffix 'Serial' is safe to be called by all processes. It immediately returns when called by the work processes. Only the master process executes it.
         * \param[in] pDataLayer init the compute layer by the metadata of the dataLayer.By the way add this dataLayer as a source to solve the compute layer.
         * \param[in] neighborFile path of the neighborhood file (text file).
         * \param[in] comptGrain 1 computeLayer cell contains comptGrain^2 dataLayer cells.
         */
        bool initSerial( RasterLayer<elemType>* pDataLayer, const char* neighborFile,int comptGrain=1);

        /**
         * \brief Init this compute layer. It should be called before solving the compute layer.
         *  
         * Metadata is also initialised by the first element of _pDataLayers
         *  Methods with suffix 'Serial' is safe to be called by all processes. It immediately returns when called by the work processes. Only the master process executes it.
         * \param[in] neighborFile path of the neighborhood file (text file).
         * \param[in] comptGrain 1 computeLayer cell contains comptGrain^2 dataLayer cells.
         */ 
        bool initSerial( const char* neighborFile,int comptGrain=1);

        void cleanDataLayers();
        vector<RasterLayer<elemType> *> *dataLayers();
        const vector<RasterLayer<elemType> *> *dataLayers() const;
        void addRasterLayer( RasterLayer<elemType> *dataLayer);
        void addRasterLayers(vector<RasterLayer<elemType> *> dataLayers);

        /**
         * \brief Get the workBR of this process.
         * \param[in] dcmpType Decomposition strategy.
         * \param[in] nSubSpcs How many parcels to divide into.
         * \param[out] subWorkBR result BR of this process.
         */ 
        bool getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR );

        /**
         * \brief Read an existing compute load file, then the compute layer can be decomposed right away.
         * \param[in] loadFile input path of an existing compute load file.
         */ 
        bool readComputeLoadFile( const char *loadFile);

        /**
         * \brief Write the compute load file out. It requires a solved compute layer.
         * \param[in] outputfile output file path.
         */         
        bool writeComputeIntensityFileSerial( const char *outputfile );
    public:
        vector< RasterLayer<elemType> *> _vDataLayers; ///< trans from this data domain to compute domain
    private:
        double _comptGrain; ///< granularity
    };
};

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer()
    :RasterLayer<elemType>() {}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer( const string RasterLayerName )
    :RasterLayer<elemType>( RasterLayerName ) {}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer( RasterLayer<elemType> * dataLayer, const int tmpGrain, const string RasterLayerName )
	: RasterLayer<elemType>( RasterLayerName ), _comptGrain( tmpGrain ) {
		_vDataLayers.push_back(dataLayer);

}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
ComputeLayer( vector<RasterLayer<elemType> *> dataLayers, const int tmpGrain, const string RasterLayerName )
    : RasterLayer<elemType>( RasterLayerName ),
      _vDataLayers( dataLayers ), _comptGrain( tmpGrain ) {
}

template<class elemType>
inline GPRO::ComputeLayer<elemType>::
~ComputeLayer() {
    //Destructor of base class is called implicitly. Only need to delete _pDataLayers.
    cleanDataLayers();
}

template<class elemType>
void GPRO::ComputeLayer<elemType>::
cleanDataLayers() {
    //dataLayers store the pointers to the layers. Deleting this vector does not influence the layers to which is pointed.
    //These layers persist for further use, until the computation finishes and they are destructed.
    vector<RasterLayer<elemType> *> vTemp;
    vTemp.swap( _vDataLayers );
}

template<class elemType>
inline vector<GPRO::RasterLayer<elemType> *> *GPRO::ComputeLayer<elemType>::
dataLayers() {
    return &_vDataLayers;
}

template<class elemType>
inline void GPRO::ComputeLayer<elemType>::
addRasterLayer(RasterLayer<elemType> *dataLayer) {
    _vDataLayers.push_back(dataLayer);
}

template<class elemType>
inline void GPRO::ComputeLayer<elemType>::
addRasterLayers(vector<RasterLayer<elemType> *> dataLayers) {
    _vDataLayers=dataLayers;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
initSerial(const char* neighborFile,int comptGrain) {
    // It is a SERIAL function. Only invoked by process 0.
    // Implicitly using members from base class is valid in Visual Studio but not allowed in gc++. i.e. _pMetaData = new MetaData() arises an error.
    
    if(GetRank()!=0) {
        return true;
    }
    if ( _vDataLayers.empty() ) {
        return false;
    }
    
    RasterLayer<elemType>::readNeighborhood(neighborFile);

    const MetaData &rhs = *( _vDataLayers[0]->_pMetaData );

    RasterLayer<elemType>::_pMetaData = new MetaData();
    MetaData *&pMetaData = RasterLayer<elemType>::_pMetaData; //Pointer as a reference. No need to delete/free.

    pMetaData->cellSize = rhs.cellSize * comptGrain;
    pMetaData->row = rhs._localworkBR.nRows() / comptGrain;
    pMetaData->row += ( rhs._localworkBR.nRows() % comptGrain ) ? 1 : 0;
    pMetaData->column = rhs._localworkBR.nCols() / comptGrain;
    pMetaData->column += ( rhs._localworkBR.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->format = rhs.format;
    pMetaData->projection = rhs.projection;
    pMetaData->noData = rhs.noData;
    pMetaData->myrank = rhs.myrank;
    pMetaData->processor_number = GetRank();
    pMetaData->_domDcmpType = rhs._domDcmpType; //ought to be NON_DCMP
    SpaceDims sdim( pMetaData->row, pMetaData->column );
    pMetaData->_glbDims = sdim;
    if ( pMetaData->_domDcmpType == NON_DCMP ) {
        CoordBR _glbWorkBR;
        RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims );
        pMetaData->_localworkBR = _glbWorkBR;
        CellCoord nwCorner( 0, 0 );
        CellCoord seCorner( pMetaData->_glbDims.nRows() - 1, pMetaData->_glbDims.nCols() - 1 );
        CoordBR subMBR( nwCorner, seCorner );
        pMetaData->_MBR = subMBR;
        pMetaData->_localdims = pMetaData->_glbDims;
    } else {
        cerr << "not support computeLayer parallel construct now." << endl;
        return false; 
    }

    pMetaData->dataType = RasterLayer<elemType>::getGDALType();

    for ( int i = 0; i < 6; i++ ) {
        pMetaData->pTransform[i] = rhs.pTransform[i];
    }
    pMetaData->pTransform[0] += rhs._localworkBR.minICol() * rhs.cellSize;
    pMetaData->pTransform[3] -= rhs._localworkBR.minIRow() * rhs.cellSize;
    pMetaData->pTransform[1] *= comptGrain;//one-pixel distance in the we/ns direction, need update
    pMetaData->pTransform[5] *= comptGrain;

    RasterLayer<elemType>::newCellSpace( pMetaData->_localdims, 0 ); //allocate

    return true;
}
template<class elemType>
bool GPRO::ComputeLayer<elemType>::
initSerial(RasterLayer<elemType>* pDataLayer, const char* neighborFile,int comptGrain) {
    // It is a SERIAL function. Only invoked by process 0.
    // Implicitly using members from base class is valid in Visual Studio but not allowed in gc++. i.e. _pMetaData = new MetaData() arises an error.

    if(GetRank()!=0) {
        return true;
    }

    if ( pDataLayer == nullptr ) {
        return false;
    }
    addRasterLayer(pDataLayer);
    RasterLayer<elemType>::readNeighborhood(neighborFile);

    MetaData &rhs = *( pDataLayer->_pMetaData );

    RasterLayer<elemType>::_pMetaData = new MetaData();
    MetaData *&pMetaData = RasterLayer<elemType>::_pMetaData; //Pointer as a reference. No need to delete/free.

    _comptGrain=comptGrain;
    pMetaData->cellSize = rhs.cellSize * comptGrain;
    //pMetaData->row = rhs._localworkBR.nRows() / comptGrain;
    //pMetaData->row += ( rhs._localworkBR.nRows() % comptGrain ) ? 1 : 0;
    //pMetaData->column = rhs._localworkBR.nCols() / comptGrain;
    //pMetaData->column += ( rhs._localworkBR.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->row = rhs._glbDims.nRows() / comptGrain;
    pMetaData->row += ( rhs._glbDims.nRows() % comptGrain ) ? 1 : 0;
    pMetaData->column = rhs._glbDims.nCols() / comptGrain;
    pMetaData->column += ( rhs._glbDims.nCols() % comptGrain ) ? 1 : 0;
    pMetaData->format = rhs.format;
    pMetaData->projection = rhs.projection;
    pMetaData->noData = rhs.noData;
    pMetaData->myrank = rhs.myrank;
    pMetaData->processor_number = GetRank(); 
    pMetaData->_domDcmpType = NON_DCMP; 
    SpaceDims sdim( pMetaData->row, pMetaData->column );
    pMetaData->_glbDims = sdim;
    if ( pMetaData->_domDcmpType == NON_DCMP ) {
        CoordBR _glbWorkBR;
        RasterLayer<elemType>::_pNbrhood->calcWorkBR( _glbWorkBR, pMetaData->_glbDims );
        pMetaData->_localworkBR = _glbWorkBR;
        CellCoord nwCorner( 0, 0 );
        CellCoord seCorner( pMetaData->_glbDims.nRows() - 1, pMetaData->_glbDims.nCols() - 1 );
        CoordBR subMBR( nwCorner, seCorner );
        pMetaData->_MBR = subMBR;
        pMetaData->_localdims = pMetaData->_glbDims;
    } else {
        cerr << "not support computeLayer parallel construct now." << endl;
        return false; 
    }

    pMetaData->dataType = RasterLayer<elemType>::getGDALType();

    for ( int i = 0; i < 6; i++ ) {
        pMetaData->pTransform[i] = rhs.pTransform[i];
    }
    pMetaData->pTransform[0] += rhs._localworkBR.minICol() * rhs.cellSize;
    pMetaData->pTransform[3] -= rhs._localworkBR.minIRow() * rhs.cellSize;
    pMetaData->pTransform[1] *= comptGrain;
    pMetaData->pTransform[5] *= comptGrain;

    RasterLayer<elemType>::newCellSpace( pMetaData->_localdims, 0 );

    return true;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
getCompuLoad( DomDcmpType dcmpType, const int nSubSpcs, CoordBR &subWorkBR ) {
    int myRank, process_nums;
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &process_nums );

    int *pDcmpIdx = new int[nSubSpcs * 4];    //the decomposition result, to be broadcast
    if ( myRank == 0 ) {    //serial, done by the master process, later broadcast to work processes
        vector<CoordBR> vDcmpBR;    //nSubSpcs parcels; MPI doesn't support custom data type, so it's useless
        
        //it should have been invisible
        vector<CoordBR> vComptDcmpBR;
        cout << "hold " << RasterLayer<elemType>::_pMetaData->_glbDims << endl;
        DeComposition<elemType> deComp( RasterLayer<elemType>::_pMetaData->_glbDims, *( RasterLayer<elemType>::_pNbrhood ));
        if ( dcmpType == ROWWISE_DCMP ) {
            deComp.valRowDcmp( vComptDcmpBR, *this, nSubSpcs );    //divide this layer by its intensity values. output to vComptDcmpBR.
        } else {
            cerr << "computLayer L388: not support until now." << endl;
        }
        //map the result to virtual BR of data domain
        CoordBR glbWorkBR;
        const Neighborhood<elemType> *pDataNbrhood = _vDataLayers[0]->nbrhood();
        pDataNbrhood->calcWorkBR( glbWorkBR, _vDataLayers[0]->_pMetaData->_glbDims );    //global workBR of data domain
        int subBegin = glbWorkBR.minIRow(), subEnd = glbWorkBR.minIRow() - 1;
        int outputRows=_vDataLayers[0]->_pMetaData->_glbDims.nRows();
        int loadFileRows=this->metaData()->_glbDims.nRows();
        _comptGrain=outputRows/double(loadFileRows);//wyj 2019-12-6 this line is added to self-adapt grain
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
        //�洢ת��������
        CellCoord nwCorner( subEnd + 1, glbWorkBR.minICol());
        CellCoord seCorner( glbWorkBR.maxIRow(), glbWorkBR.maxICol());
        CoordBR subMBR( nwCorner, seCorner );
        vDcmpBR.push_back( subMBR );
        pDcmpIdx[4 * i] = subEnd + 1;
        pDcmpIdx[4 * i + 1] = glbWorkBR.minICol();
        pDcmpIdx[4 * i + 2] = glbWorkBR.maxIRow();
        pDcmpIdx[4 * i + 3] = glbWorkBR.maxICol();
    }

    MPI_Barrier( MPI_COMM_WORLD );    //work processes wait here to get their BRs.
    MPI_Bcast( pDcmpIdx, process_nums * 4, MPI_INT, 0, MPI_COMM_WORLD );
    CellCoord nwCorner2( pDcmpIdx[4 * myRank], pDcmpIdx[4 * myRank + 1] );
    CellCoord seCorner2( pDcmpIdx[4 * myRank + 2], pDcmpIdx[4 * myRank + 3] );
    CoordBR tmpWorkBR( nwCorner2, seCorner2 );
    subWorkBR = tmpWorkBR;    //this line is tested, OK.
    //cout<<"tmpWorkBR "<<tmpWorkBR<<" subWorkBR "<<subWorkBR<<" shoule be the same"<<endl;
    MPI_Barrier( MPI_COMM_WORLD );
    delete[]pDcmpIdx;

    return true;
}
template<class elemType>
bool GPRO::ComputeLayer<elemType>::
readComputeLoadFile(const char *loadFile) {
    if(GetRank()!=0){
        //MPI_Barrier( MPI_COMM_WORLD );
        return true;
    }
    
    RasterLayer<elemType>::readFile(loadFile,NON_DCMP);
    return true;
}

template<class elemType>
bool GPRO::ComputeLayer<elemType>::
writeComputeIntensityFileSerial( const char *outputfile ) {
    GDALAllRegister();
    if(GetRank()!=0) {
        MPI_Barrier( MPI_COMM_WORLD );
        return true;
    }
    if ( !RasterLayer<elemType>::createFile( outputfile )) {
        cout << "create file is not correct!" << endl;
        MPI_Finalize();
    }

    GDALDataset *poDataset = NULL;
    poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
    if ( poDataset == NULL) {
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