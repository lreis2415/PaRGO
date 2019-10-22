#ifndef RasterLayer_H
#define RasterLayer_H

/***************************************************************************
* RasterLayer.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::RasterLayer
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
#include "deComposition.h"
#include <fstream>
#include <sstream>

#include"mpi.h"
#include<gdal_priv.h>
using namespace std;

namespace GPRO
{
  template <class elemType>
  class RasterLayer
  {
    public:
      RasterLayer();
      RasterLayer(const string RasterLayerName = "Untitled");
      RasterLayer(const RasterLayer<elemType> &rhs);
      ~RasterLayer();

      RasterLayer<elemType>& operator=(const RasterLayer<elemType> &rhs);

      const string &name() const;
      void name(const string &RasterLayerName);
     
      int id() const;
      const string title() const;
     

      void cleanCellSpace();
      void cleanNbrhood();
      void cleanMetaData();

      bool hasCellSpace() const;
      bool hasNbrhood() const;
	  bool hasMetaData() const;

     
      bool newCellSpace();
      bool newCellSpace(const SpaceDims& dims);
      bool newCellSpace(const SpaceDims& dims,
                        const elemType &initVal);
      bool newCellSpace(int nRows, int nCols);
      bool newCellSpace(int nRows, int nCols,
                        const elemType &initVal);

      bool newNbrhood();
      bool newNbrhood(const vector<CellCoord> &vNbrCoords,
                      double weight = 1.0);
      bool newNbrhood(const vector<CellCoord> &vNbrCoords,
                      const vector<double> &vNbrWeights);
                      
      template<class elemType2>
      bool newNbrhood(const Neighborhood<elemType2> &nbr);

     

      CellSpace<elemType> *cellSpace();
      const CellSpace<elemType> *cellSpace() const;
      Neighborhood<elemType> *nbrhood();
      const Neighborhood<elemType> *nbrhood() const;
      MetaData *metaData();

	  bool copyLayerInfo(const RasterLayer<elemType> &rhs);

	  GDALDataType getType();
	  MPI_Datatype getMPIType();
	  //IO function
	  bool readNeighborhood(const char* neighborfile);
	  bool readFile(const char* inputfile);
	  bool writeFile(const char* outputfile);
	  bool createFile(const char *outputfile);


  public:
	  MetaData* _pMetaData;

    protected:
      string _name;
      CellSpace<elemType> *_pCellSpace; // 
      Neighborhood<elemType> *_pNbrhood; // 
	  
      
  };
};

template <class elemType>
inline GPRO::RasterLayer<elemType>::
RasterLayer()
  :
   _name("Untitled"),
   _pCellSpace(0),
   _pNbrhood(0)
   {}

template <class elemType>
inline GPRO::RasterLayer<elemType>::
RasterLayer(const string RasterLayerName)
  :_name(RasterLayerName),
   _pCellSpace(0),
   _pNbrhood(0)
   //_pMetaData(0)
{
  
}



template <class elemType>
inline GPRO::RasterLayer<elemType>::
~RasterLayer() 
{
  cleanCellSpace();
  cleanNbrhood();
  cleanMetaData();
}

template <class elemType>
inline GPRO::RasterLayer<elemType>& GPRO::RasterLayer<elemType>::
operator=(const RasterLayer<elemType> &rhs) 
{
  if(this == &rhs)
  {
    return *this;
  }
  
  /*
  if(_name == "Untitled") {
    _name = rhs._name + "_copy";
  }
  */
  if(rhs._pCellSpace) 
  {
    if(_pCellSpace)
	{
      *(_pCellSpace) = *(rhs._pCellSpace);
    }
    else 
	{
      _pCellSpace = new CellSpace<elemType>(*(rhs._pCellSpace));
    }
  }
  else 
  {
    cleanCellSpace();
  }
  
  if(rhs._pNbrhood)
  {
    if(_pNbrhood)
	{
      *(_pNbrhood) = *(rhs._pNbrhood);
    }
    else 
	{
      _pNbrhood = new Neighborhood<elemType>(*(rhs._pNbrhood));
    }
  }
  else 
  {
    cleanNbrhood();
  }

  return *this;
}

template <class elemType>
inline const string& GPRO::RasterLayer<elemType>::
name() const
{
  return _name;
}

template <class elemType>
inline void GPRO::RasterLayer<elemType>::
name(const string &RasterLayerName) 
{
  _name = RasterLayerName;
}



template <class elemType>
inline int GPRO::RasterLayer<elemType>::
id() const 
{
  int myID = -1;
  
  return myID;
}

template <class elemType>
inline const string GPRO::RasterLayer<elemType>::
title() const 
{
  ostringstream myTitle;
  myTitle << _name << id();
  return myTitle.str();
}



template <class elemType>
void GPRO::RasterLayer<elemType>::
cleanCellSpace() 
{
  if(_pCellSpace) 
  {
    delete _pCellSpace;
    _pCellSpace = 0;
  }
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
cleanNbrhood() 
{
  if(_pNbrhood) 
  {
    delete _pNbrhood;
    _pNbrhood = 0;
  }
}

template <class elemType>
void GPRO::RasterLayer<elemType>::
cleanMetaData() 
{
  if(_pMetaData) 
  {
    delete _pMetaData;
    _pMetaData = 0;
  }
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
hasCellSpace() const
{
  bool hasIt = true;
  if(!_pCellSpace ||
     _pCellSpace->empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no CellSpace associated with RasterLayer[" \
         << title() << "]" << endl;
    hasIt = false;
  }
  return hasIt;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
hasNbrhood() const 
{
  bool hasIt = true;
  if(!_pNbrhood) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: no neighborhood associated with RasterLayer[" \
         << title() << "]" << endl;
    hasIt = false;
  }
  return hasIt;
}



template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace() 
{
 
    cleanCellSpace();
  
    _pCellSpace = new CellSpace<elemType>();
    if(!_pCellSpace) 
	{
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to new a CellSpace" \
           << endl;
      return false;
    }
  
  return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(const SpaceDims& dims) 
{
  
	cleanCellSpace();
 
	_pCellSpace = new CellSpace<elemType>(dims);
	if(!_pCellSpace) 
	{
		cerr << __FILE__ << " " << __FUNCTION__ \
			<< " Error: unable to new a CellSpace" \
			<< endl;
		return false;
	}
  
  return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(const SpaceDims& dims,
             const elemType &initVal) 
{
 
    cleanCellSpace();

    _pCellSpace = new CellSpace<elemType>(dims, initVal);
    if(!_pCellSpace)
	{
		cerr << __FILE__ << " " << __FUNCTION__ \
			<< " Error: unable to new a CellSpace" \
			<< endl;
		return false;
    }
  
  return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(int nRows, int nCols)
{
  return newCellSpace(SpaceDims(nRows, nCols));
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newCellSpace(int nRows, int nCols,
             const elemType &initVal)
{
  return newCellSpace(SpaceDims(nRows, nCols), initVal);
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood() 
{
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>();
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood(const vector<CellCoord> &vNbrCoords,
           double weight) 
{
 
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>(vNbrCoords, weight);
  if(!_pNbrhood) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}

template <class elemType>
inline GPRO::CellSpace<elemType>* GPRO::RasterLayer<elemType>::
cellSpace() 
{
  return _pCellSpace;
}

template <class elemType>
inline const GPRO::CellSpace<elemType>* GPRO::RasterLayer<elemType>::
cellSpace() const 
{
  return _pCellSpace;
}

template <class elemType>
inline GPRO::Neighborhood<elemType>* GPRO::RasterLayer<elemType>::
nbrhood() 
{
  return _pNbrhood;
}

template <class elemType>
inline const GPRO::Neighborhood<elemType>* GPRO::RasterLayer<elemType>::
nbrhood() const 
{
  return _pNbrhood;
}


template <class elemType>
bool GPRO::RasterLayer<elemType>::
newNbrhood(const vector<CellCoord> &vNbrCoords,
           const vector<double> &vNbrWeights) 
{
  
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>(vNbrCoords, vNbrWeights);
  if(!_pNbrhood) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}

template<class elemType> template<class elemType2>
bool GPRO::RasterLayer<elemType>::
newNbrhood(const Neighborhood<elemType2> &nbr) 
{
  
  cleanNbrhood();
  _pNbrhood = new Neighborhood<elemType>(nbr);
  if(!_pNbrhood) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to new a Neighborhood on process[" \
         << title() << "]" << endl;
    return false;
  }
  return true;
}
   

template <class elemType>
MPI_Datatype  GPRO::RasterLayer<elemType>::
getMPIType()
{
	elemType style;
	if( typeid(style) == typeid(float) )
	{
		return MPI_FLOAT;
	}
	else if ( typeid(style) == typeid(double) )
	{
		return MPI_DOUBLE;
		//cout<<"double"<<endl;
	}
	else if ( typeid(style) == typeid(int) )
	{
		return MPI_INT;
		//cout<<"int"<<endl;
	}
	else if ( typeid(style) == typeid(unsigned int) )
	{
		return MPI_UNSIGNED;
		//cout<<"unsigned int"<<endl;
	}
	else if ( typeid(style) == typeid(char) )
	{
		return MPI_CHAR;
		//cout<<"char"<<endl;
	}
	else
	{
		return MPI_BYTE;
	}
}

template <class elemType>
GDALDataType GPRO::RasterLayer<elemType>::
getType()
{
	elemType style;
	GDALDataType dataType;
	//cout<<typeid(style).name()<<endl;
	if( typeid(style) == typeid(float) )
	{
		dataType = GDT_Float32;
		//cout<<"float"<<endl;
	}
	else if ( typeid(style) == typeid(double) )
	{
		dataType = GDT_Float64;
		//cout<<"double"<<endl;
	}
	else if ( typeid(style) == typeid(int) )
	{
		dataType = GDT_Int32;
		//cout<<"int"<<endl;
	}
	else if ( typeid(style) == typeid(unsigned int) )
	{
		dataType = GDT_UInt32;
		//cout<<"unsigned int"<<endl;
	}
	else if ( typeid(style) == typeid(char) )
	{
		dataType = GDT_Byte;
		//cout<<"char"<<endl;
	}
	else
	{
		dataType = GDT_Unknown;
	}
	return dataType;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
readNeighborhood(const char* neighborfile)
{
	newNbrhood();
	
	fstream nbrFile(neighborfile, ios::in);
	nbrFile >> (*_pNbrhood); /* Initialize the Neighborhood object by loading the neighbors stored in a ASCII file */
	nbrFile.close();
	return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
copyLayerInfo(const RasterLayer<elemType> &rhs)
{
	_pMetaData = new MetaData();
	if(_pMetaData == NULL)
	{
		//do something
		cout<<"[ERROR] MetaData is not allocate correct"<<endl;
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
	_pMetaData->dataType = getType();

	for(int i = 0; i < 6; i++)
	{
		_pMetaData->pTransform[i] = rhs._pMetaData->pTransform[i];
	}
	
	newNbrhood(*(rhs.nbrhood())); //allocate and init
	newCellSpace(_pMetaData->_localdims,_pMetaData->noData); //allocate
	//newCellSpace(_pMetaData->_localdims,0); //allocate

	return true;
}

template <class elemType>
bool GPRO::RasterLayer<elemType>::
readFile(const char* inputfile)
{
	GDALAllRegister();
	//cout<<inputfile<<endl;
	GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(inputfile, GA_ReadOnly );
	//GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(inputfile, GA_ReadOnly );
	if( poDatasetsrc == NULL /*检查是否正常打开文件*/)
	{
		//do something
		cout<<"[ERROR] data file is not open correct"<<endl;
		exit(1);
	}
	_pMetaData = new MetaData();

	if(_pMetaData == NULL)
	{
		//do something
		cout<<"[ERROR] MetaData is not allocate correct"<<endl;
		exit(1);
	}

	
	poDatasetsrc->GetGeoTransform(_pMetaData->pTransform);
	_pMetaData->projection = poDatasetsrc->GetProjectionRef();
	GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand( 1 );

	_pMetaData->noData = poBandsrc->GetNoDataValue();
	_pMetaData->row = poBandsrc->GetYSize();
	_pMetaData->column = poBandsrc->GetXSize();
	SpaceDims sdim(_pMetaData->row, _pMetaData->column);
	_pMetaData->_glbDims = sdim;
	_pMetaData->cellSize = _pMetaData->pTransform[1];
	_pMetaData->format = "GTiff";

	MPI_Comm_rank(MPI_COMM_WORLD, &_pMetaData->myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &_pMetaData->processor_number);


	DeComposition<elemType> deComp(_pMetaData->_glbDims, *_pNbrhood);
	deComp.rowDcmp(*_pMetaData, _pMetaData->processor_number);
	//SpaceDims dims(_metadata.subRow, _metadata.column);
	newCellSpace(_pMetaData->_localdims);
	//_pCellSpace
	///*if (_matrix == NULL)
	//{
	//	cout<<"[ERROR] the allocation of memory is error!"<<endl;
	//	MPI_Finalize();
	//	return false;
	//}*/
	
	//cout<<"My rank is "<<myrank<<"! temp is "<<temp<<" column is "<<column<<" subRow is "<<subRow<<endl;
	_pMetaData->dataType = getType();

	//cout<<"_pMetaData->myrank "<<_pMetaData->myrank<<"  _pMetaData->row "<<_pMetaData->row<<endl;
	//cout<<"_pMetaData->myrank "<<_pMetaData->myrank<<"  _pMetaData->_MBR.minIRow() "<<_pMetaData->_MBR.minIRow()<<"  _pMetaData->_localdims.nRows()  "<<_pMetaData->_localdims.nRows()<<"  _pNbrhood->minIRow()  "<<_pNbrhood->minIRow()<<endl;
	
	poBandsrc->RasterIO(GF_Read, 0, _pMetaData->_MBR.minIRow(), _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pCellSpace->_matrix, _pMetaData->_localdims.nCols(), _pMetaData->_localdims.nRows(), _pMetaData->dataType, 0, 0);

	if (poDatasetsrc != NULL)
	{
		//
		//poDatasetsrc->~GDALDataset();
		GDALClose((GDALDatasetH)poDatasetsrc);
		poDatasetsrc = NULL;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return true;
}


template <class elemType>
bool GPRO::RasterLayer<elemType>::
createFile(const char* outputfile)
{
	GDALAllRegister();
	if(_pMetaData->myrank == 0)
	{
		GDALDriver* poDriver = NULL;
		poDriver = GetGDALDriverManager()->GetDriverByName(_pMetaData->format.c_str());
		if (poDriver == NULL)
		{
			cout<<"poDriver is NULL."<<endl;
			return false;
		}
		char **papszMetadata = NULL;
		papszMetadata = CSLSetNameValue(papszMetadata, "BLOCKXSIZE", "256");
		papszMetadata = CSLSetNameValue( papszMetadata, "BLOCKYSIZE", "1" );

		if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) )
			;
			//cout<< "This format driver supports Create() method."<<endl;
		if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) )
			;
			//cout<< "This format driver supports CreateCopy() method."<<endl;
	
		GDALDataset* poDataset = poDriver->Create(outputfile, _pMetaData->_glbDims.nCols(),  _pMetaData->_glbDims.nRows(), 1, _pMetaData->dataType, papszMetadata);
		if (poDataset == NULL)
		{
			//do something
			cout<<"poDatasetdest is NULL"<<endl;
			return false;
		}
		poDataset->SetGeoTransform( _pMetaData->pTransform );	
		poDataset->SetProjection(_pMetaData->projection.c_str());

		if (poDataset != NULL)
		{
			GDALClose((GDALDatasetH)poDataset);
			poDataset = NULL;
		}
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return true;
}


template <class elemType>
bool GPRO::RasterLayer<elemType>::
writeFile(const char* outputfile)
{
	GDALAllRegister();
	
	if(!createFile(outputfile))
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
	if(_pMetaData->myrank == 0)
	{
		poBanddest->SetNoDataValue(_pMetaData->noData);
	}

	if(_pMetaData->processor_number == 1)
	{
		//cout<<"_pMetaData->myrank "<<_pMetaData->myrank <<"  _pMetaData->_glbDims.nCols() "<<_pMetaData->_glbDims.nCols()<<"   _pMetaData->_glbDims.nRows() "<<_pMetaData->_glbDims.nRows()<<endl;
			
		poBanddest->RasterIO(GF_Write, 0, 0, _pMetaData->_glbDims.nCols(),_pMetaData->_glbDims.nRows(), _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nRows(), _pMetaData->dataType, 0, 0);
	}
	else
	{
		if(_pMetaData->myrank == 0)
		{
			//cout<<"_pMetaData->myrank "<<_pMetaData->myrank <<"  _pMetaData->_localworkBR.maxIRow()+1 "<<_pMetaData->_localworkBR.maxIRow()+1<<endl;
		
			poBanddest->RasterIO(GF_Write, 0, 0, _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+1, _pCellSpace->_matrix, _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+1, _pMetaData->dataType, 0, 0);
		}
		else if(_pMetaData->myrank == ( _pMetaData->processor_number - 1) )
		{
			//poBanddest->RasterIO(GF_Write, 0,  _pMetaData->_MBR.minIRow()-_pNbrhood->minIRow(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+1, _pCellSpace->_matrix-_pNbrhood->minIRow(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+1, _pMetaData->dataType, 0, 0);
			//cout<<"_pMetaData->myrank "<<_pMetaData->myrank <<"  _pMetaData->_MBR.minIRow()-_pNbrhood->minIRow() "<<_pMetaData->_MBR.minIRow()-_pNbrhood->minIRow()<<"  _pMetaData->_localworkBR.maxIRow()+1 "<<_pMetaData->_localworkBR.maxIRow()+1<<endl;
			poBanddest->RasterIO(GF_Write, 0,  _pMetaData->_MBR.minIRow()-_pNbrhood->minIRow(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+1, _pCellSpace->_matrix-_pNbrhood->minIRow()*_pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+1, _pMetaData->dataType, 0, 0);
			
		}
		else
		{
			poBanddest->RasterIO(GF_Write, 0,  _pMetaData->_MBR.minIRow()-_pNbrhood->minIRow(), _pMetaData->_glbDims.nCols(), _pMetaData->_localworkBR.maxIRow()+_pNbrhood->minIRow()+1, _pCellSpace->_matrix-_pNbrhood->minIRow()*_pMetaData->_glbDims.nCols(), _pMetaData->_glbDims.nCols(),  _pMetaData->_localworkBR.maxIRow()+_pNbrhood->minIRow()+1, _pMetaData->dataType, 0, 0);
		}
	}
	
	if (poDataset != NULL)
	{
		GDALClose((GDALDatasetH)poDataset);
		poDataset = NULL;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return true;
}



#endif
