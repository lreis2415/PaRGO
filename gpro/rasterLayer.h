/**
 * \file basicTypes
 * \author Qin ChengZhi (Qincz@lreis.ac.cn), Wang Yujing
 * \brief Header file for class GPRO::RasterLayer
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


#ifndef RasterLayer_H
#define RasterLayer_H

#include "basicTypes.h"
#include "cellSpace.h"
#include "neighborhood.h"
#include "metaData.h"
#include "deComposition.h"
#include "utility.h"
#include <fstream>
#include <sstream>

#include "mpi.h"
#include <gdal_priv.h>

using namespace std;
#define EPS 0.0000001

namespace GPRO {
    /**
     * \ingroup gpro
     * \class RasterLayer
     * \brief To store a raster layer.
     */
    template <class elemType>
    class RasterLayer {
    public:

        /**
         *constructors
         */
        RasterLayer();
        RasterLayer(string layerName);
        RasterLayer(const RasterLayer<elemType>& rhs);
        virtual ~RasterLayer();

        /**
         *copy constructor
         */
        RasterLayer<elemType>& operator=(const RasterLayer<elemType>& rhs);

        /**
         *getters and setters
         */
        const string& name() const;
        void name(const string& layerName);
        unsigned int id() const;
        void id(int layerID);
        const string title() const;

        /**
         *getters
         */
        CellSpace<elemType>* cellSpace();
        const CellSpace<elemType>* cellSpace() const;
        Neighborhood<elemType>* nbrhood();
        const Neighborhood<elemType>* nbrhood() const;
        MetaData* metaData();
        GDALDataType getGDALType();
        MPI_Datatype getMPIType();

        /**
         *Destruction functions
         */
        void cleanCellSpace();
        void cleanNbrhood();
        void cleanMetaData();

        bool hasCellSpace() const;
        bool hasNbrhood() const;
        bool hasMetaData() const;

        /**
         *init CellSpace
         */
        bool newCellSpace();
        bool newCellSpace(const SpaceDims& dims);
        bool newCellSpace(const SpaceDims& dims, const elemType& initVal);
        bool newCellSpace(int nRows, int nCols);
        bool newCellSpace(int nRows, int nCols, const elemType& initVal);

        /**
         *init Neighborhood
         */
        bool newNbrhood();
        bool newLocalNbrhood();
        bool newNbrhood(const vector<CellCoord>& vNbrCoords, double weight = 1.0);
        bool newRectangleNbrhood(int radiusExcludingCenter);
        bool newNbrhood(const vector<CellCoord>& vNbrCoords, const vector<double>& vNbrWeights);
        template <class elemType2>
        bool newNbrhood(const Neighborhood<elemType2>& nbr);

        /**
         * \brief Initialize the Neighborhood object by loading the neighbors stored in a ASCII file
         * \param[in] ASCII file
         */
        bool readNeighborhood(const char* neighborfile);
        /**
         * \brief Initialize the Neighborhood object by loading the neighbors stored in a ASCII file
         * \param[in] ASCII file
         */
        bool readNeighborhoodSerial(const char* neighborfile);
        /**
         * \brief read raster layer.
         *  
         * Include sync. Should be called by all processes.
         *  First calculate the desired area according to dcmpType; then each process reads the appointed area. 
         * \param[in] inputfile raster file
         * \param[in] dcmpType decomposition type
         */
        bool readFile(const char* inputfile, DomDcmpType dcmpType = NON_DCMP);

        bool newUpscaleFile(const char* inputfile, int g, DomDcmpType dcmpType = NON_DCMP);

        /**
         * \brief read raster layer.
         *  
         * Include sync. Should be called by all processes.
         *  First calculate the desired area according to dcmpType; then each process reads the appointed area. 
         * \param[in] inputfile raster file
         * \param[in] dcmpType decomposition type. It is not used in this method, but may be for further use.
         */
        bool readFile(const char* inputfile, const CoordBR& subWorkBR, DomDcmpType dcmpType = NON_DCMP);


        /**
         *all processes read raster layer according to the decomposition type. 
         */
        // bool allReadFile( const char *inputfile, DomDcmpType dcmpType = NON_DCMP );


        /**
         * \brief read the whole raster layer by the master process.
         *  
         * Methods with suffix 'Serial' is safe to be called by all processes. It immediately return when called by the work processes.
         * \param[in] inputfile raster file
         */
        bool readGlobalFileSerial(const char* inputfile);

        /**
         * \brief read the information and init _pMetaData
         * \param[in] poDatasetsrc data source
         * \param[in] dcmpType decomposition type. info should be coordinate with data.
         */
        GDALRasterBand* readFileInfo(GDALDataset* poDatasetsrc, DomDcmpType dcmpType);

        GDALRasterBand* readUpscaleFileInfo(GDALDataset* poDatasetsrc, DomDcmpType dcmpType, int g);

        /**
         * \brief decompose this layer and init _pMetaData & CellSpace.
         * \param[in] dcmpType decomposition type.
         */
        bool layerDcmp(DomDcmpType dcmpType);

        /**
         * \brief decompose this layer and init _pMetaData & CellSpace.
         * \param[in] subWorkBR by which the members are initialised.
         */
        bool layerDcmp(const CoordBR& subWorkBR);

        /**
         * \brief copy the metadata and cellspace of input.
         * \param[in] rhs The input RasterLayer. Abbr for 'right hand side'.
         */
        bool copyLayerInfo(const RasterLayer<elemType>& rhs);

        /**
         * \brief copy the metadata of input.
         * \param[in] rhs The input RasterLayer. Abbr for 'right hand side'.
         */
        bool copyLayerMetadata(const RasterLayer<elemType>& rhs);

        //void col2row( int &subRows, int &lastSubRows );	//invert the col-wise data to row-wise

        /**
         * \brief output this RasterLayer to raster file
         *  
         * according to the DCMP_TYPE in MetaData.
         *  support Row and Col output now
         * \param[in] rhs The input RasterLayer.
         */
        bool writeFile(const char* outputfile);


        int rowAtOtherLayer(RasterLayer<double>* layer, int row);
        int rowAtOtherLayer(RasterLayer<int>* layer, int row);
        int colAtOtherLayer(RasterLayer<double>* layer, int col);
        int colAtOtherLayer(RasterLayer<int>* layer, int col);
        bool isNodata(elemType value);
    protected:

        /**
         * \brief create a file.
         *   
         * It must be called by all processes, but the file will only be created by the main process serially
         * \param[in] outputfile The input RasterLayer.
         */
        bool createFile(const char* outputfile);

        /**
         * \brief parallel output, row-wisely
         * \param[in] outputfile The input RasterLayer.
         */
        bool rowWriteFile(const char* outputfile);

        /**
         * \brief parallel IO, refer to Qin13_TransactionsInGIS
         * \param[in] outputfile The input RasterLayer.
         */
        bool colWriteFile(const char* outputfile);

    public:
        MetaData* _pMetaData; ///< raster layer metadata

    protected:
        string _strLayerName; ///< layer name
        CellSpace<elemType>* _pCellSpace; ///< raster data stored in this process
        Neighborhood<elemType>* _pNbrhood; ///< neighbor data
        int _nLayerID; ///< layer ID
    };
};

template <class elemType>
GPRO::RasterLayer<elemType>::
RasterLayer()
    :
    _pMetaData(nullptr),
    _strLayerName("Untitled"),
    _pCellSpace(nullptr),
    _pNbrhood(nullptr),
    _nLayerID(-1) {
    _pNbrhood = new Neighborhood<elemType>(true);
}

template <class elemType>
GPRO::RasterLayer<elemType>::
RasterLayer(const string layerName)
    : _pMetaData(nullptr),
      _strLayerName(layerName),
      _pCellSpace(nullptr),
      _pNbrhood(nullptr),
      _nLayerID(-1) {
    _pNbrhood = new Neighborhood<elemType>(true);
}

template <class elemType>
GPRO::RasterLayer<elemType>::
~RasterLayer() {
    cleanCellSpace();
    cleanNbrhood();
    cleanMetaData();
}

template <class elemType>
GPRO::RasterLayer<elemType>& GPRO::RasterLayer<elemType>::
operator=(const RasterLayer<elemType>& rhs) {
    if (this == &rhs) {
        return *this;
    }

    if (!_strLayerName.empty()) {
        _strLayerName = rhs._strLayerName + "_copy";
    }

    ////Todo: copy LayerID

    if (rhs._pCellSpace) {
        if (_pCellSpace) {
            cleanCellSpace();
            *(_pCellSpace) = *(rhs._pCellSpace);
        }
        else {
            _pCellSpace = new CellSpace<elemType>(*(rhs._pCellSpace));
        }
    }
    else {
        cleanCellSpace();
    }

    if (rhs._pNbrhood) {
        if (_pNbrhood) {
            cleanNbrhood();
            *(_pNbrhood) = *(rhs._pNbrhood);
        }
        else {
            _pNbrhood = new Neighborhood<elemType>(*(rhs._pNbrhood));
        }
    }
    else {
        cleanNbrhood();
    }

    return *this;
}

template <class elemType>
const string& GPRO::RasterLayer<elemType>::
name() const {
    return _strLayerName;
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
name(const string& layerName) {
    _strLayerName = layerName;
}

template <class elemType>
unsigned int GPRO::RasterLayer<elemType>::
id() const {
    return _nLayerID;
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
id(const int layerID) {
    if (layerID < 0) {
        cout << "The layerID should be bigger than zero." << endl;
    }
    _nLayerID = layerID;
}

template <class elemType>
const string GPRO::RasterLayer<elemType>::
title() const {
    ostringstream myTitle;
    return myTitle.str();
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
isNodata(elemType value) {
    return value - _pMetaData->noData < EPS;
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
cleanCellSpace() {
    if (_pCellSpace) {
        delete _pCellSpace;
        _pCellSpace = 0;
    }
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
cleanNbrhood() {
    if (_pNbrhood) {
        delete _pNbrhood;
        _pNbrhood = 0;
    }
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
cleanMetaData() {
    if (_pMetaData) {
        delete _pMetaData;
        _pMetaData = nullptr;
    }
}


template <class elemType>
bool GPRO::RasterLayer<elemType>::
hasCellSpace() const {
    bool hasIt = true;
    if (!_pCellSpace || _pCellSpace->empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: no CellSpace associated with RasterLayer["
            << title() << "]" << endl;
        hasIt = false;
    }

    return hasIt;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
hasNbrhood() const {
    bool hasIt = true;
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: no neighborhood associated with RasterLayer["
            << title() << "]" << endl;
        hasIt = false;
    }

    return hasIt;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
hasMetaData() const {
    bool hasIt = true;
    if (!_pMetaData) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: no metaData associated with RasterLayer["
            << title() << "]" << endl;
        hasIt = false;
    }

    return hasIt;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace() {
    cleanCellSpace();

    _pCellSpace = new CellSpace<elemType>();
    if (!_pCellSpace) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a CellSpace"
            << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(const SpaceDims& dims) {
    cleanCellSpace();
    _pCellSpace = new CellSpace<elemType>(dims);
    if (!_pCellSpace) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a CellSpace"
            << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(const SpaceDims& dims, const elemType& initVal) {
    cleanCellSpace();
    _pCellSpace = new CellSpace<elemType>(dims, initVal);
    if (!_pCellSpace) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a CellSpace"
            << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(int nRows, int nCols) {
    return newCellSpace(SpaceDims(nRows, nCols));
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(int nRows, int nCols, const elemType& initVal) {
    return newCellSpace(SpaceDims(nRows, nCols), initVal);
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood() {
    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>();
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a Neighborhood on process["
            << title() << "]" << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newLocalNbrhood() {
    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>(true);
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a Neighborhood on process["
            << title() << "]" << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood(const vector<CellCoord>& vNbrCoords, double weight) {

    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>(vNbrCoords, weight);
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a Neighborhood on process["
            << title() << "]" << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newRectangleNbrhood(int radiusExcludingCenter) {
    int r = radiusExcludingCenter;
    cleanNbrhood();
    vector<CellCoord> nbr;
    for (int i = -r; i <= r; ++i) {
        for (int j = -r; j <= r; ++j) {
            CellCoord c(i, j);
            nbr.push_back(c);
        }
    }
    _pNbrhood = new Neighborhood<elemType>(nbr, 1);
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a Neighborhood on process["
            << title() << "]" << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood(const vector<CellCoord>& vNbrCoords, const vector<double>& vNbrWeights) {
    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>(vNbrCoords, vNbrWeights);
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a Neighborhood on process["
            << title() << "]" << endl;
        return false;
    }

    return true;
}

template <class elemType>
template <class elemType2>
bool GPRO::RasterLayer<elemType>::
newNbrhood(const Neighborhood<elemType2>& nbr) {
    cleanNbrhood();
    _pNbrhood = new Neighborhood<elemType>(nbr);
    if (!_pNbrhood) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to new a Neighborhood on process["
            << title() << "]" << endl;
        return false;
    }

    return true;
}

template <class elemType>
GPRO::CellSpace<elemType>* GPRO::RasterLayer<elemType>::
cellSpace() {
    return _pCellSpace;
}

template <class elemType>
const GPRO::CellSpace<elemType>* GPRO::RasterLayer<elemType>::
cellSpace() const {
    return _pCellSpace;
}

template <class elemType>
GPRO::Neighborhood<elemType>* GPRO::RasterLayer<elemType>::
nbrhood() {
    return _pNbrhood;
}

template <class elemType>
const GPRO::Neighborhood<elemType>* GPRO::RasterLayer<elemType>::
nbrhood() const {
    return _pNbrhood;
}

template <class elemType>
GPRO::MetaData* GPRO::RasterLayer<elemType>::
metaData() {
    return _pMetaData;
}

template <class elemType>
MPI_Datatype GPRO::RasterLayer<elemType>::
getMPIType() {
    elemType style;
    if (typeid(style) == typeid(float)) {
        return MPI_FLOAT;
    }
    if (typeid(style) == typeid(double)) {
        return MPI_DOUBLE;
    }
    if (typeid(style) == typeid(int)) {
        return MPI_INT;
    }
    if (typeid(style) == typeid(unsigned int)) {
        return MPI_UNSIGNED;
    }
    if (typeid(style) == typeid(char)) {
        return MPI_CHAR;
    }
    else {
        return MPI_BYTE;
    }
}

template <class elemType>
GDALDataType GPRO::RasterLayer<elemType>::
getGDALType() {
    elemType style;
    GDALDataType dataType;
    if (typeid(style) == typeid(float)) {
        dataType = GDT_Float32;
    }
    else if (typeid(style) == typeid(double)) {
        dataType = GDT_Float64;
    }
    else if (typeid(style) == typeid(int)) {
        dataType = GDT_Int32;
    }
    else if (typeid(style) == typeid(unsigned int)) {
        dataType = GDT_UInt32;
    }
    else if (typeid(style) == typeid(char)) {
        dataType = GDT_Byte;
    }
    else {
        dataType = GDT_Unknown;
    }
    return dataType;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
readNeighborhood(const char* neighborfile) {
    newNbrhood();
    fstream nbrFile(neighborfile, ios::in);
    nbrFile >> (*_pNbrhood);
    nbrFile.close();

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
readNeighborhoodSerial(const char* neighborfile) {
    if (GetRank() != 0) {
        return true;
    }
    newNbrhood();
    fstream nbrFile(neighborfile, ios::in);
    nbrFile >> (*_pNbrhood);
    nbrFile.close();

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
readGlobalFileSerial(const char* inputfile) {
    if (GetRank() != 0) {
        //MPI_Barrier( MPI_COMM_WORLD );
        return true;
    }
    GDALAllRegister();
    GDALDataset* poDatasetsrc = static_cast<GDALDataset*>(GDALOpen(inputfile, GA_ReadOnly));

    GDALRasterBand* poBandsrc = readFileInfo(poDatasetsrc, NON_DCMP);

    layerDcmp(NON_DCMP);

    poBandsrc->RasterIO(GF_Read, 0, _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pMetaData->dataType, 0, 0);

    if (poDatasetsrc != nullptr) {
        GDALClose(static_cast<GDALDatasetH>(poDatasetsrc));
        poDatasetsrc = nullptr;
    }
    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
copyLayerMetadata(const RasterLayer<elemType>& rhs) {
    _pMetaData = new MetaData();
    if (_pMetaData == nullptr) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit(1);
    }

    _pMetaData->cellSize = rhs._pMetaData->cellSize;
    _pMetaData->row = rhs._pMetaData->row;
    _pMetaData->column = rhs._pMetaData->column;
    _pMetaData->format = rhs._pMetaData->format;
    _pMetaData->projection = rhs._pMetaData->projection;
    _pMetaData->noData = rhs._pMetaData->noData;
    _pMetaData->myrank = rhs._pMetaData->myrank;
    _pMetaData->processor_number = rhs._pMetaData->processor_number;
    _pMetaData->_domDcmpType = rhs._pMetaData->_domDcmpType;
    _pMetaData->_glbDims = rhs._pMetaData->_glbDims;
    _pMetaData->_MBR = rhs._pMetaData->_MBR;
    _pMetaData->_localdims = rhs._pMetaData->_localdims;
    _pMetaData->_localworkBR = rhs._pMetaData->_localworkBR;
    _pMetaData->dataType = getGDALType();

    for (int i = 0; i < 6; i++) {
        _pMetaData->pTransform[i] = rhs._pMetaData->pTransform[i];
    }
    newNbrhood(*(rhs.nbrhood())); //allocate and init
    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
copyLayerInfo(const RasterLayer<elemType>& rhs) {
    copyLayerMetadata(rhs);
    newCellSpace(_pMetaData->_localdims, rhs._pMetaData->noData);
    return true;
}

//TODO: overload copyLayerInfo to copy from a different elemType.

template <class elemType>
GDALRasterBand* GPRO::RasterLayer<elemType>::
readFileInfo(GDALDataset* poDatasetsrc, DomDcmpType dcmpType) {
    if (poDatasetsrc == nullptr) {
        cout << "[ERROR] data file is not open correct" << endl;
        exit(1);
    }

    _pMetaData = new MetaData();
    if (_pMetaData == nullptr) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit(1);
    }

    poDatasetsrc->GetGeoTransform(_pMetaData->pTransform);
    _pMetaData->projection = poDatasetsrc->GetProjectionRef();
    GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand(1);

    _pMetaData->noData = poBandsrc->GetNoDataValue();
    _pMetaData->row = poBandsrc->GetYSize();
    _pMetaData->column = poBandsrc->GetXSize();
    _pMetaData->cellSize = _pMetaData->pTransform[1];
    _pMetaData->format = "GTiff";
    _pMetaData->_domDcmpType = dcmpType;
    _pMetaData->dataType = getGDALType();
    SpaceDims sdim(_pMetaData->row, _pMetaData->column);
    _pMetaData->_glbDims = sdim;

    return poBandsrc;
}

template <class elemType>
GDALRasterBand* GPRO::RasterLayer<elemType>::
readUpscaleFileInfo(GDALDataset* poDatasetsrc, DomDcmpType dcmpType, int g) {
    if (poDatasetsrc == nullptr) {
        cout << "[ERROR] data file is not open correct" << endl;
        exit(1);
    }

    _pMetaData = new MetaData();
    if (_pMetaData == nullptr) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit(1);
    }

    poDatasetsrc->GetGeoTransform(_pMetaData->pTransform);
    _pMetaData->projection = poDatasetsrc->GetProjectionRef();
    GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand(1);

    _pMetaData->noData = poBandsrc->GetNoDataValue();
    _pMetaData->row = poBandsrc->GetYSize() / g;
    _pMetaData->column = poBandsrc->GetXSize() / g;

    _pMetaData->pTransform[1] *= g; //w-e pixel resolution
    _pMetaData->pTransform[5] *= g; //n-s pixel resolution
    _pMetaData->cellSize = _pMetaData->pTransform[1];

    _pMetaData->format = "GTiff";
    _pMetaData->_domDcmpType = dcmpType;
    _pMetaData->dataType = getGDALType();
    SpaceDims sdim(_pMetaData->row, _pMetaData->column);
    _pMetaData->_glbDims = sdim;

    return poBandsrc;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
layerDcmp(DomDcmpType dcmpType) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_pMetaData->myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &_pMetaData->processor_number);

    DeComposition<elemType> deComp(_pMetaData->_glbDims, *_pNbrhood);
    if (ROWWISE_DCMP == _pMetaData->_domDcmpType) {
        deComp.rowDcmp(*_pMetaData, _pMetaData->processor_number);
    }
    else {
        if (COLWISE_DCMP == _pMetaData->_domDcmpType) {
            deComp.colDcmp(*_pMetaData, _pMetaData->processor_number);
        }
        else {
            if (BLOCK_DCMP == _pMetaData->_domDcmpType) {
                cout << __FILE__ << " " << __FUNCTION__
                    << "Error: not support this dcmpType_" << dcmpType
                    << " right now" << endl; //TODO
                return false;
            }
            else {
                if (NON_DCMP == _pMetaData->_domDcmpType) {
                    _pMetaData->_localdims = _pMetaData->_glbDims;
                    CoordBR _glbWorkBR;
                    _pNbrhood->calcWorkBR(_glbWorkBR, _pMetaData->_glbDims);
                    _pMetaData->_localworkBR = _glbWorkBR;
                    CellCoord nwCorner(0, 0);
                    CellCoord seCorner(_pMetaData->_glbDims.nRows() - 1, _pMetaData->_glbDims.nCols() - 1);
                    CoordBR subMBR(nwCorner, seCorner);
                    _pMetaData->_MBR = subMBR;
                }
                else {
                    cout << __FILE__ << " " << __FUNCTION__
                        << "Error: not support this dcmpType_" << dcmpType
                        << " right now" << endl;
                    return false;
                }
            }
        }
    }
    newCellSpace(_pMetaData->_localdims);
    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
layerDcmp(const CoordBR& subWorkBR) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_pMetaData->myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &_pMetaData->processor_number);

    int glbBegin = subWorkBR.nwCorner().iRow();
    int glbEnd = subWorkBR.seCorner().iRow();
    glbEnd = min(glbEnd, _pMetaData->_glbDims.nRows()-1);
    CellCoord nwCorner(glbBegin + _pNbrhood->minIRow(), 0);
    CellCoord seCorner(glbEnd + _pNbrhood->maxIRow(), _pMetaData->_glbDims.nCols() - 1);
    CoordBR subMBR(nwCorner, seCorner);
    _pMetaData->_MBR = subMBR;
    SpaceDims dims(subMBR.nRows(), subMBR.nCols());
    _pMetaData->_localdims = dims;

    //subWorkBR is a global index, while _localworkBR is a local index
    CoordBR workBR;
    if (!_pNbrhood->calcWorkBR(workBR, dims)) {
        return false;
    }
    _pMetaData->_localworkBR = workBR;

    newCellSpace(_pMetaData->_localdims);
    return true;
}

//Todo: test irregular neighbor
template <class elemType>
bool GPRO::RasterLayer<elemType>::
readFile(const char* inputfile, DomDcmpType dcmpType) {
    GDALAllRegister();
    GDALDataset* poDatasetsrc = static_cast<GDALDataset*>(GDALOpen(inputfile, GA_ReadOnly));

    GDALRasterBand* poBandsrc = readFileInfo(poDatasetsrc, dcmpType);

    if (!layerDcmp(dcmpType)) return false;

    poBandsrc->RasterIO(GF_Read, _pMetaData->_MBR.minICol(), _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                        _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                        _pMetaData->dataType, 0, 0);

    if (poDatasetsrc != nullptr) {
        GDALClose(static_cast<GDALDatasetH>(poDatasetsrc));
        poDatasetsrc = nullptr;
    }
    // MPI_Barrier( MPI_COMM_WORLD ); // wyj: maybe useless.

    return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newUpscaleFile(const char* inputfile, int g, DomDcmpType dcmpType) {
    GDALAllRegister();
    GDALDataset* poDatasetsrc = static_cast<GDALDataset*>(GDALOpen(inputfile, GA_ReadOnly));

    GDALRasterBand* poBandsrc = readUpscaleFileInfo(poDatasetsrc, dcmpType, g);

    if (!layerDcmp(dcmpType)) return false;
    newCellSpace(_pMetaData->_localdims, _pMetaData->noData);

    if (poDatasetsrc != nullptr) {
        GDALClose(static_cast<GDALDatasetH>(poDatasetsrc));
        poDatasetsrc = nullptr;
    }

    return true;
}

//Todo: CoordBR should be invisible to user, may use RowID or ColID instead.
//is dcmpType still in use?
template <class elemType>
bool GPRO::RasterLayer<elemType>::
readFile(const char* inputfile, const CoordBR& subWorkBR, DomDcmpType dcmpType) {
    //why cleanMetaData crashes£¿
    if (_pMetaData) {
        cleanMetaData();
    } // to overwrite instead of clean is ok ?
    GDALAllRegister();
    GDALDataset* poDatasetsrc = static_cast<GDALDataset*>(GDALOpen(inputfile, GA_ReadOnly));
    GDALRasterBand* poBandsrc = readFileInfo(poDatasetsrc, dcmpType);

    if (!layerDcmp(subWorkBR)) return false;

    poBandsrc->RasterIO(GF_Read, 0, _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                        _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(),
                        _pMetaData->dataType, 0, 0);

    if (poDatasetsrc != nullptr) {
        GDALClose(static_cast<GDALDatasetH>(poDatasetsrc));
        poDatasetsrc = nullptr;
    }
    MPI_Barrier(MPI_COMM_WORLD); // if needed?

    return true;
}

//template<class elemType>
//void GPRO::RasterLayer<elemType>::
//col2row( int &subRows, int &lastSubRows ){}

// master process creates file
template <class elemType>
bool GPRO::RasterLayer<elemType>::
createFile(const char* outputfile) {
    GDALAllRegister();

    if (_pMetaData->myrank == 0) {
        GDALDriver* poDriver = nullptr;
        poDriver = GetGDALDriverManager()->GetDriverByName(_pMetaData->format.c_str());
        if (poDriver == nullptr) {
            cout << "poDriver is NULL." << endl;
            return false;
        }
        char** papszMetadata = nullptr;
        // papszMetadata = CSLSetNameValue(papszMetadata, "BLOCKXSIZE", "256");
        // papszMetadata = CSLSetNameValue(papszMetadata, "BLOCKYSIZE", "1");
        // papszMetadata = CSLSetNameValue(papszMetadata, "COMPRESS", "LZW");

        if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE));
        //cout<< "This format driver supports Create() method."<<endl;
        if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATECOPY, FALSE));
        //cout<< "This format driver supports CreateCopy() method."<<endl;

        GDALDataset* poDataset = poDriver->Create(outputfile, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                                                  1, _pMetaData->dataType, papszMetadata);
        if (poDataset == nullptr) {
            cout << "poDatasetdest is NULL" << endl;
            return false;
        }
        poDataset->SetGeoTransform(_pMetaData->pTransform);
        poDataset->SetProjection(_pMetaData->projection.c_str());

        if (poDataset != nullptr) {
            GDALClose(static_cast<GDALDatasetH>(poDataset));
            poDataset = nullptr;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return true;
}


//TODO: test irregular neighbor, like -1/0/1/2
template <class elemType>
bool GPRO::RasterLayer<elemType>::
writeFile(const char* outputfile) {
    if (COLWISE_DCMP == _pMetaData->_domDcmpType) {
        return colWriteFile(outputfile);
    }
    else if (BLOCK_DCMP == _pMetaData->_domDcmpType) {
        cout << __FILE__ << " " << __FUNCTION__
            << "Error: not support this dcmpType_" << _pMetaData->_domDcmpType
            << " right now" << endl; //unfinished
        return false;
    }
    else if (NON_DCMP == _pMetaData->_domDcmpType) {
        //done by master or not to support?
    }
    else {
        return rowWriteFile(outputfile);
    }
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
rowWriteFile(const char* outputfile) {
    GDALAllRegister();

    if (!createFile(outputfile)) {
        cout << __FILE__ << " " << __FUNCTION__
            << " Error: create file failed!" << endl;
        MPI_Finalize();
    }
    GDALDataset* poDataset = static_cast<GDALDataset*>(GDALOpen(outputfile, GA_Update));
    if (poDataset == nullptr) {
        cout << "data file is not open correct" << endl;
        exit(1);
    }

    GDALRasterBand* poBanddest = poDataset->GetRasterBand(1);
    if (poBanddest == nullptr) {
        cout << "poBanddest is NULL" << endl;
        exit(1);
    }
    if (_pMetaData->myrank == 0) {
        poBanddest->SetNoDataValue(_pMetaData->noData); //why only master process does this
    }
    if (_pMetaData->processor_number == 1) {
        poBanddest->RasterIO(GF_Write, 0, 0, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                             _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                             _pMetaData->dataType, 0, 0);
    }
    else {
        int nXOff=0;
        int nYOff;
        int nXSize=_pMetaData->_glbDims.nCols();
        int nYSize;
        int nBufXSize=_pMetaData->_glbDims.nCols();
        int nBufYSize;
        elemType* pData;
        if (_pMetaData->myrank == 0) {
            nYOff=0;
            nYSize=_pMetaData->_localworkBR.maxIRow() + 1;
            pData=_pCellSpace->_matrix;
            nBufYSize=_pMetaData->_localworkBR.maxIRow() + 1;
        }
        else if (_pMetaData->myrank == (_pMetaData->processor_number - 1)) {
            nYOff=_pMetaData->_MBR.minIRow() - _pNbrhood->minIRow();
            nYSize=_pMetaData->_localworkBR.maxIRow() + 1;
            int mirow= _pNbrhood->minIRow();
            pData=_pCellSpace->_matrix - _pNbrhood->minIRow() * _pMetaData->_glbDims.nCols();
            nBufYSize=_pMetaData->_localworkBR.maxIRow() + 1;
        }
        else {
            nYOff=_pMetaData->_MBR.minIRow() - _pNbrhood->minIRow();
            nYSize= _pMetaData->_localworkBR.maxIRow() + _pNbrhood->minIRow() + 1;
            pData=_pCellSpace->_matrix - _pNbrhood->minIRow() * _pMetaData->_glbDims.nCols();
            nBufYSize=_pMetaData->_localworkBR.maxIRow() + _pNbrhood->minIRow() + 1;
        }
        cout<<*_pCellSpace;
        poBanddest->RasterIO(GF_Write, nXOff, nYOff, nXSize, nYSize, pData, nBufXSize, nBufYSize, _pMetaData->dataType, 0, 0);
    }

    MPI_Barrier(MPI_COMM_WORLD); //if needed?
    GDALClose(poDataset);
    poDataset = nullptr;
    // if (poDataset != nullptr) {
    //     GDALClose(static_cast<GDALDatasetH>(poDataset));
    //     poDataset = nullptr;
    // }
    return true;
}


template <class elemType>
bool GPRO::RasterLayer<elemType>::
colWriteFile(const char* outputfile) {
    int myRank, process_nums;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
    MPI_Datatype datatype = getMPIType();

    GDALAllRegister();
    if (!createFile(outputfile)) {
        cout << "create file is not correct!" << endl;
        MPI_Finalize();
    }

    GDALDataset* poDataset = nullptr;
    poDataset = static_cast<GDALDataset*>(GDALOpen(outputfile, GA_Update));
    if (poDataset == nullptr) {
        cout << "data file is not open correct" << endl;
        exit(1);
    }
    GDALRasterBand* poBanddest = poDataset->GetRasterBand(1);
    if (poBanddest == nullptr) {
        cout << "poBanddest is NULL" << endl;
        exit(1);
    }
    if (_pMetaData->myrank == 0) {
        poBanddest->SetNoDataValue(_pMetaData->noData);
    }

    if (_pMetaData->processor_number == 1) {
        poBanddest->RasterIO(GF_Write, 0, 0, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                             _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(),
                             _pMetaData->dataType, 0, 0);
    }
    else {
        //phase 1, re-distribute
        int subRows = 0, lastSubRows = 0;
        if (_pMetaData->_glbDims.nRows() % process_nums == 0) {
            subRows = _pMetaData->_glbDims.nRows() / process_nums;
            lastSubRows = subRows;
        }
        else {
            subRows = _pMetaData->_glbDims.nRows() / process_nums + 1;
            lastSubRows = _pMetaData->_glbDims.nRows() - (_pMetaData->_glbDims.nRows() / process_nums + 1) * (process_nums - 1);
        }
        int glbWorkCols = _pMetaData->_glbDims.nCols() + _pNbrhood->minICol() - _pNbrhood->maxICol();
        // allocate memory, store itself, and store right after receiving one
        int cellCount = 0;
        if (myRank < process_nums - 1) {
            cellCount = subRows * _pMetaData->_glbDims.nCols();
        }
        else {
            cellCount = lastSubRows * _pMetaData->_glbDims.nCols();
        }

        elemType* pRowMatrix = new elemType[cellCount];
        for (int rRow = 0; rRow < subRows; ++rRow) {
            if (myRank == 0) {
                for (int rCol = 0; rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol) {
                    pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                }
            }
            else {
                if (myRank == process_nums - 1) {
                    for (int rCol = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                         rCol <= _pMetaData->_MBR.maxICol() + _pNbrhood->maxICol(); ++rCol) {
                        pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] =
                            _pCellSpace->_matrix[(rRow + myRank * subRows) * _pMetaData->_localdims.nCols() + rCol - _pMetaData->_MBR.minICol()];
                    }
                }
                else {
                    for (int rCol = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                         rCol <= _pMetaData->_MBR.maxICol(); ++rCol) {
                        pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] =
                            _pCellSpace->_matrix[(rRow + myRank * subRows) * _pMetaData->_localdims.nCols() + rCol - _pMetaData->_MBR.minICol()];
                    }
                }
            }
        }
        if (nullptr == pRowMatrix) {
            cout << "c " << myRank << endl;
        }
        MPI_Status status;
        //The size of rankID's myRank-th parcel is cols*RowOfMyRank, create/init/update/delete this memmory
        //Q Isn't the size fixed? Why communicate the size of temp space?
        for (int objRank = 0; objRank < process_nums; ++objRank) {
            if (objRank == myRank) {
                continue;
            }
            //communicate objRank's workBR.nCols() and MBR.minCol() to myRank
            int sendnCols, sendColPos;
            if (myRank == 0) {
                sendnCols = _pMetaData->_localworkBR.nCols() + 1;
                sendColPos = 0;
            }
            else {
                if (myRank == process_nums - 1) {
                    sendnCols = _pMetaData->_localworkBR.nCols() + 1;
                    sendColPos = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                }
                else {
                    sendnCols = _pMetaData->_localworkBR.nCols();
                    sendColPos = _pMetaData->_MBR.minICol() - _pNbrhood->minICol();
                }
            }
            int recvnCols, recvColPos;
            MPI_Send(&sendnCols, 1, MPI_INT, objRank, 1, MPI_COMM_WORLD);
            MPI_Recv(&recvnCols, 1, MPI_INT, objRank, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&sendColPos, 1, MPI_INT, objRank, 2, MPI_COMM_WORLD);
            MPI_Recv(&recvColPos, 1, MPI_INT, objRank, 2, MPI_COMM_WORLD, &status);
            int sendMatrixSize, recvMatrixSize;
            if (myRank < process_nums - 1) {
                recvMatrixSize = recvnCols * subRows;
            }
            else {
                recvMatrixSize = recvnCols * lastSubRows;
            }
            if (objRank < process_nums - 1) {
                sendMatrixSize = sendnCols * subRows;
            }
            else {
                sendMatrixSize = sendnCols * lastSubRows;
            }
            elemType* sendMatrix = new elemType[sendMatrixSize];
            elemType* recvMatrix = new elemType[recvMatrixSize];
            //myRank will receive from rankID the size of "my localRows * rankID's workBR rows" (get by communication)
            //the data to send and receive should be stored in a temporary array, then communicate and map
            cellCount = 0;
            for (int rRow = objRank * subRows; rRow < (objRank + 1) * subRows; ++rRow) {
                if (rRow >= _pMetaData->_glbDims.nRows()) {
                    break;
                }
                if (myRank == 0) {
                    for (int rCol = 0; rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol) {
                        sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                    }
                }
                else {
                    if (myRank == process_nums - 1) {
                        for (int rCol = _pMetaData->_localworkBR.minICol(); rCol <= _pMetaData->_localworkBR.maxICol() + _pNbrhood->maxICol(); ++rCol) {
                            sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                        }
                    }
                    else {
                        for (int rCol = _pMetaData->_localworkBR.minICol(); rCol <= _pMetaData->_localworkBR.maxICol(); ++rCol) {
                            sendMatrix[cellCount++] = _pCellSpace->_matrix[rRow * _pMetaData->_localdims.nCols() + rCol];
                        }
                    }
                }
            }

            int styleSize;
            MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &styleSize);
            int bufsize = MPI_BSEND_OVERHEAD + styleSize * (sendMatrixSize + recvMatrixSize);
            void* buf = malloc(bufsize);
            MPI_Buffer_attach(buf, bufsize);

            MPI_Bsend(sendMatrix, sendMatrixSize, datatype, objRank, 1, MPI_COMM_WORLD);
            MPI_Recv(recvMatrix, recvMatrixSize, datatype, objRank, 1, MPI_COMM_WORLD, &status);
            cellCount = 0;
            //write the recvMatrix into a matrix
            for (int rRow = 0; rRow < subRows; ++rRow) {
                if (myRank == process_nums - 1 && rRow >= lastSubRows) {
                    break;
                }
                for (int rCol = recvColPos; rCol < recvColPos + recvnCols; ++rCol) {
                    pRowMatrix[rRow * _pMetaData->_glbDims.nCols() + rCol] = recvMatrix[cellCount++];
                }
            }
            //if( myRank == 1 )
            //cout<<recvColPos<<" "<<recvCols<<endl;
            MPI_Buffer_detach(&buf, &bufsize);
            MPI_Buffer_attach(buf, bufsize); //not need?
        }

        if (myRank < process_nums - 1) {
            poBanddest->RasterIO(GF_Write, 0, myRank * subRows, _pMetaData->_glbDims.nCols(), subRows,
                                 pRowMatrix, _pMetaData->_glbDims.nCols(), subRows,
                                 _pMetaData->dataType, 0, 0);
        }
        else {
            cout << "a " << myRank << " " << pRowMatrix << endl;
            poBanddest->RasterIO(GF_Write, 0, myRank * subRows, _pMetaData->_glbDims.nCols(), lastSubRows,
                                 pRowMatrix, _pMetaData->_glbDims.nCols(), lastSubRows,
                                 _pMetaData->dataType, 0, 0);
        }
    }
    cout << "b " << myRank << endl;

    MPI_Barrier(MPI_COMM_WORLD); //will wait forever if not executed by all processes
    if (poDataset != nullptr) {
        GDALClose(static_cast<GDALDatasetH>(poDataset));
        poDataset = nullptr;
    }

    return true;
}

template <class elemType>
int GPRO::RasterLayer<elemType>::
rowAtOtherLayer(RasterLayer<double>* layer, int row) {
    double scale = static_cast<double>(layer->metaData()->_glbDims.nRows()) / _pMetaData->_glbDims.nRows();
    return scale * row;
}

template <class elemType>
int GPRO::RasterLayer<elemType>::
rowAtOtherLayer(RasterLayer<int>* layer, int row) {
    double scale = static_cast<double>(layer->metaData()->_glbDims.nRows()) / _pMetaData->_glbDims.nRows();
    return scale * row;
}

template <class elemType>
int GPRO::RasterLayer<elemType>::
colAtOtherLayer(RasterLayer<double>* layer, int col) {
    double scale = static_cast<double>(layer->metaData()->_glbDims.nCols()) / _pMetaData->_glbDims.nCols();
    return scale * col;
}

template <class elemType>
int GPRO::RasterLayer<elemType>::
colAtOtherLayer(RasterLayer<int>* layer, int col) {
    double scale = static_cast<double>(layer->metaData()->_glbDims.nCols()) / _pMetaData->_glbDims.nCols();
    return scale * col;
}

#endif
