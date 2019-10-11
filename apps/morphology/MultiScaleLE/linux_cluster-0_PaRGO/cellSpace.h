#ifndef CELLSPACE_H
#define CELLSPACE_H

/***************************************************************************
* cellSpace.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::Transition, and GPRO::CellSpace
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
#include <string>
#include <iostream>
#include <typeinfo>
#include"mpi.h"
#include<gdal_priv.h>


using namespace std;

namespace GPRO 
{
  template <class elemType>
  class CellSpace;

  template <class elemType>
  class Neighborhood;

  template <class elemType>
  class Transition 
  {
    public:
      Transition(bool onlyUpdtCtrCell = true,
                 bool needFinalize = true,
                 bool needExchange = true,
                 bool edgesFirst = true)
        :_pCellSpace(0),
         _pNbrhood(0),
         _onlyUpdtCtrCell(onlyUpdtCtrCell),
         _needFinalize(needFinalize),
         _needExchange(needExchange),
         _edgesFirst(edgesFirst) {}

      virtual ~Transition() {}

      bool onlyUpdtCtrCell() const 
	  {
        return _onlyUpdtCtrCell;
      }
      bool needFinalize() const 
	  {
        return _needFinalize;
      }
      bool needExchange() const 
	  {
        return _needExchange; 
      }
      bool edgesFirst() const 
	  {
        return _edgesFirst;
      }
      void needFinalize(bool ndFnlz) 
	  {
        _needFinalize = ndFnlz;
      }

      virtual bool cellSpace(CellSpace<elemType> *pCellSpc) 
	  {
        _pCellSpace = pCellSpc;
        if(!_pCellSpace) 
		{
          return false;
        }
        return true;
      }
      virtual bool nbrhood(Neighborhood<elemType> *pNbrhood) 
	  {
        _pNbrhood = pNbrhood;
        return true;
      }

      virtual bool evaluate(vector<pair<int, elemType> > &vUpdtedCells,
                            const CellCoord &coord) {}

      virtual bool finalize(const elemType &val,
                            const CellCoord &coord) {}

      virtual int workload(const CoordBR &workBR) {}

    protected:
      CellSpace<elemType> *_pCellSpace;
      Neighborhood<elemType> *_pNbrhood;
      bool _onlyUpdtCtrCell;
      bool _needFinalize;
      bool _needExchange;
      bool _edgesFirst;
  };

  template <class elemType>
  ostream& operator<<(ostream &os,
                      const CellSpace<elemType> &cellSpace);

  template <class elemType>
  istream& operator>>(istream &is,
                      CellSpace<elemType> &cellSpace);

  template <class elemType>
  class CellSpace
  {
    public:
      // Constructor
      CellSpace();
      CellSpace(const SpaceDims& dims);
      CellSpace(const SpaceDims& dims,
                const elemType &initVal);
      CellSpace(int nRows, int nCols);
      CellSpace(int nRows, int nCols, 
                const elemType &initVal);
      CellSpace(const CellSpace<elemType> &rhs);

      // Destructor
      ~CellSpace();
	  //read write create file!
	  //////////////////////////
	 /* bool readFile(const char* inputfile);
	  bool writeFile(const char* outputfile);
	  bool createFile(const char *outputfile);*/

      bool initMem(const SpaceDims& dims);
      bool initVals(const elemType &initVal);
      void clear();
      
      template <class elemType2>
      bool equalDim(const CellSpace<elemType2> &rhs) const;
      bool empty() const;
      int nRows() const;
      int nCols() const;
      int size() const;
      const SpaceDims& dims() const; 
      bool validCoord(const CellCoord& coord,
                      bool warning = true) const;
      bool validCoord(int iRow, int iCol,
                      bool warning = true) const;
      bool validIdx(int idx,
                    bool warning = true) const;

      elemType* operator[](int iRow);
      const elemType* operator[](int iRow) const;
      CellSpace<elemType>& operator=(const CellSpace<elemType> &rhs);

      int coord2idx(const CellCoord& coord) const;
      int coord2idx(int iRow, int iCol) const;
      const CellCoord idx2coord(int idx) const;

      bool values(vector<elemType> &vVals) const;

      bool find(IntVect &vFoundIdxs,
                const elemType &val) const;
      bool find(IntVect &vFoundIdxs,
                const vector<elemType> &vVals) const;
      bool find(IntVect &vFoundIdxs,
                const vector<elemType> &vVals,
                const CoordBR &rectangle) const;

      int count(const elemType &val) const;
      int count(const elemType &val,
                const CoordBR &rectangle) const;
      int count(const elemType &val,
                const IntVect &vExcldIdxs) const;
      int count(const elemType &val,
                const CoordBR &rectangle,
                const IntVect &vExcldIdxs) const;

      int count(const vector<elemType> &vVals) const;
      int count(const vector<elemType> &vVals,
                const CoordBR &rectangle) const;
      int count(const vector<elemType> &vVals,
                const IntVect &vExcldIdxs) const;
      int count(const vector<elemType> &vVals,
                const CoordBR &rectangle,
                const IntVect &vExcldIdxs) const;

      template<class Predicate>
      int count(Predicate pred) const;

      template<class Predicate>
      int count(Predicate pred,
                const CoordBR &rectangle) const;
      
      bool updateRowCol(const vector<elemType> &vNewVals,
                        int iDim,
                        int iRowCol);
      bool updateCells(const vector< BasicCell<elemType> > &vCells);
      bool updateCells(const vector<pair<int, elemType> > &vCells);
      bool update(Transition<elemType> *pTransition,
                  Neighborhood<elemType> *pNbrhood,
                  const CoordBR *const pWorkBR = 0);
      bool update(IntVect &vIdxs2Eval,
                  Transition<elemType> *pTransition,
                  Neighborhood<elemType> *pNbrhood,
                  const CoordBR *const pWorkBR = 0);
      bool updateFinalize(Transition<elemType> *pTransition = 0,
                          Neighborhood<elemType> *pNbrhood = 0);
      bool add2UpdtMap(const elemType &newVal,
                       const CellCoord &coord);
      void cleanUpdtMap();
	 // GDALDataType getType();

  public:
	  MetaData _metadata;
	  elemType *_matrix;
    protected:
      SpaceDims _dims;
      map<elemType, IntVect > _mUpdtCells;
  };

};

/****************************************
*             iostream Methods          *
*****************************************/
template <class elemType>
inline ostream& GPRO::
operator<<(ostream &os,
           const CellSpace<elemType> &cellSpace) 
{
  int nRows = cellSpace.nRows();
  int nCols = cellSpace.nCols();
  os << nRows << " " << nCols << endl;
  for(int iRow = 0; iRow < nRows; iRow++) 
  {
    for(int iCol = 0; iCol < nCols; iCol++) 
	{
      os << cellSpace[iRow][iCol] << " ";
    } // End of iCol loop
    os << endl;
  } // End of iRow loop
  return os;
}

template <class elemType>
inline istream& GPRO::
operator>>(istream &is,
           CellSpace<elemType> &cellSpace) 
{
  int nRows, nCols;
  is >> nRows;
  is >> nCols;
  SpaceDims dims(nRows, nCols);

  if(!cellSpace.empty() &&
     cellSpace.dims() != dims) 
  {
    cellSpace.clear();
    cellSpace.initMem(dims);
  }
  else if(cellSpace.empty()) 
  {
    cellSpace.initMem(dims);
  }

  for(int iRow = 0; iRow < nRows; iRow++) 
  {
    for(int iCol = 0; iCol < nCols; iCol++) 
	{
      int iElem = iRow*nCols + iCol;
      is >> cellSpace[iRow][iCol];
    } // End of iCol loop
  } // End of iRow loop
  return is;
}

/****************************************
*             Public Methods            *
*****************************************/ 

//template <class elemType>
//GDALDataType GPRO::CellSpace<elemType>::
//getType()
//{
//	elemType style;
//	GDALDataType dataType;
//	cout<<typeid(style).name()<<endl;
//	if( typeid(style) == typeid(float) )
//	{
//		dataType = GDT_Float32;
//		//cout<<"float"<<endl;
//	}
//	else if ( typeid(style) == typeid(double) )
//	{
//		dataType = GDT_Float64;
//		//cout<<"double"<<endl;
//	}
//	else if ( typeid(style) == typeid(int) )
//	{
//		dataType = GDT_Int32;
//		//cout<<"int"<<endl;
//	}
//	else if ( typeid(style) == typeid(unsigned int) )
//	{
//		dataType = GDT_UInt32;
//		//cout<<"unsigned int"<<endl;
//	}
//	else if ( typeid(style) == typeid(char) )
//	{
//		dataType = GDT_Byte;
//		//cout<<"char"<<endl;
//	}
//	else
//	{
//		dataType = GDT_Unknown;
//	}
//	return dataType;
//}
//
//template <class elemType>
//bool GPRO::CellSpace<elemType>::
//createFile(const char* outputfile)
//{
//	GDALAllRegister();
//	if(_metadata.myrank == 0)
//	{
//		GDALDriver* poDriver = NULL;
//		poDriver = GetGDALDriverManager()->GetDriverByName(_metadata.format.c_str());
//		if (poDriver == NULL)
//		{
//			cout<<"poDriver is NULL."<<endl;
//			return false;
//		}
//		char **papszMetadata = NULL;
//		papszMetadata = CSLSetNameValue(papszMetadata, "BLOCKXSIZE", "256");
//		papszMetadata = CSLSetNameValue( papszMetadata, "BLOCKYSIZE", "1" );
//
//		if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) )
//			;
//			//cout<< "This format driver supports Create() method."<<endl;
//		if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) )
//			;
//			//cout<< "This format driver supports CreateCopy() method."<<endl;
//	
//		GDALDataset* poDataset = poDriver->Create(outputfile, _metadata.column, _metadata.row, 1, _metadata.dataType, papszMetadata);
//		if (poDataset == NULL)
//		{
//			//do something
//			cout<<"poDatasetdest is NULL"<<endl;
//			return false;
//		}
//		poDataset->SetGeoTransform( _metadata.pTransform );	
//		poDataset->SetProjection(_metadata.projection.c_str());
//
//		if (poDataset != NULL)
//		{
//			GDALClose((GDALDatasetH)poDataset);
//			poDataset = NULL;
//		}
//		
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	return true;
//}
//
//template <class elemType>
//bool GPRO::CellSpace<elemType>::
//writeFile(const char* outputfile)
//{
//	GDALAllRegister();
//	
//	if(!createFile(outputfile))
//	{
//		cout<<"create file is not correct!"<<endl;
//		MPI_Finalize();
//	}
//	GDALDataset* poDataset = NULL;
//	poDataset = (GDALDataset *) GDALOpen( outputfile, GA_Update );
//	if( poDataset == NULL /*检查是否正常打开文件*/)
//	{
//		//do something
//		cout<<"data file is not open correct"<<endl;
//		exit(1);
//	}
//	GDALRasterBand*	poBanddest = poDataset->GetRasterBand(1);
//	if (poBanddest == NULL)
//	{
//		//do something
//		cout<<"poBanddest is NULL"<<endl;
//		exit(1);
//	}
//	if(_metadata.myrank == 0)
//	{
//		poBanddest->SetNoDataValue(_metadata.noData);
//	}
//	
//	poBanddest->RasterIO(GF_Write, 0,  _metadata.myrank*_metadata.interval, _metadata.column, _metadata.subRow, _matrix, _metadata.column, _metadata.subRow, _metadata.dataType, 0, 0);
//	
//	if (poDataset != NULL)
//	{
//		GDALClose((GDALDatasetH)poDataset);
//		poDataset = NULL;
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	return true;
//}

//template <class elemType>
//bool GPRO::CellSpace<elemType>::
//readFile(const char* inputfile)
//{
//	GDALAllRegister();
//	//cout<<inputfile<<endl;
//	GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(inputfile, GA_ReadOnly );
//	//GDALDataset* poDatasetsrc = (GDALDataset *) GDALOpen(inputfile, GA_ReadOnly );
//	if( poDatasetsrc == NULL /*检查是否正常打开文件*/)
//	{
//		//do something
//		cout<<"[ERROR] data file is not open correct"<<endl;
//		exit(1);
//	}
//	poDatasetsrc->GetGeoTransform(_metadata.pTransform);
//	_metadata.projection = poDatasetsrc->GetProjectionRef();
//	GDALRasterBand* poBandsrc = poDatasetsrc->GetRasterBand( 1 );
//
//	_metadata.noData = poBandsrc->GetNoDataValue();
//	_metadata.row = poBandsrc->GetYSize();
//	_metadata.column = poBandsrc->GetXSize();
//	_metadata.cellSize = _metadata.pTransform[1];
//	_metadata.format = "GTiff";
//	MPI_Comm_rank(MPI_COMM_WORLD, &_metadata.myrank);
//	MPI_Comm_size(MPI_COMM_WORLD, &_metadata.processor_number);
//	
//
//	_metadata.interval = (_metadata.row + _metadata.processor_number - 1)/_metadata.processor_number;
//
//	if ( _metadata.myrank == (_metadata.processor_number - 1) )
//	{
//		_metadata.subRow = _metadata.row - (_metadata.processor_number - 1)*_metadata.interval;
//		_metadata.block = _metadata.subRow*_metadata.column;
//		//_matrix = (elemType*)CPLMalloc(sizeof(elemType)*_metadata.block);
//		
//	}
//	else
//	{
//		_metadata.subRow = _metadata.interval;
//		_metadata.block = _metadata.subRow*_metadata.column;
//		//_matrix = (elemType*)CPLMalloc(sizeof(elemType)*_metadata.block);
//		
//		
//	}
//
//	SpaceDims dims(_metadata.subRow, _metadata.column);
//
//	if(!this->empty() &&
//	this->dims() != dims) 
//	{
//		this->clear();
//		this->initMem(dims);
//	}
//	else if(this->empty()) 
//	{
//		this->initMem(dims);
//	}
//
//
//	if (_matrix == NULL)
//	{
//		cout<<"[ERROR] the allocation of memory is error!"<<endl;
//		MPI_Finalize();
//		return false;
//	}
//	
//	//cout<<"My rank is "<<myrank<<"! temp is "<<temp<<" column is "<<column<<" subRow is "<<subRow<<endl;
//	_metadata.dataType = getType();
//	
//	poBandsrc->RasterIO(GF_Read, 0, _metadata.myrank*_metadata.interval, _metadata.column, _metadata.subRow, _matrix, _metadata.column, _metadata.subRow, _metadata.dataType, 0, 0);
//
//	if (poDatasetsrc != NULL)
//	{
//		//
//		//poDatasetsrc->~GDALDataset();
//		GDALClose((GDALDatasetH)poDatasetsrc);
//		poDatasetsrc = NULL;
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	return true;
//}

template <class elemType>
inline bool GPRO::CellSpace<elemType>::
initMem(const SpaceDims& dims) 
{
  bool done = true;
  if(!dims.valid())
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid Dimensions [" << dims \
         << "] to allocate memory for the CellSpace" \
         << endl;
    done = false;
  }
  else 
  {
    _dims = dims;
    _matrix = new elemType[dims.size()];
    if(!_matrix) 
	{
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: unable to allocate memory for the CellSpace" \
           << endl;
      done = false;	
    }
  }
  return done;
}

template <class elemType>
inline bool GPRO::CellSpace<elemType>::
initVals(const elemType &initVal) 
{
  bool done = true;
  if(empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to initialize an empty CellSpace" \
         << endl;
    done = false;
  }
  for(int iElem = 0; iElem < size(); iElem++) 
  {
    _matrix[iElem] = initVal;	
  }
  return done;
}

template <class elemType>
inline void GPRO::CellSpace<elemType>::
clear()
{
  if(_matrix) 
  {
    delete[] _matrix;
  }
}

template <class elemType>
inline GPRO::CellSpace<elemType>::
CellSpace()
  :_matrix(0) {}

template <class elemType>
inline GPRO::CellSpace<elemType>::
CellSpace(const SpaceDims& dims)
{
  if(!initMem(dims))
  {
    exit(-1); 
  }
}

template <class elemType>
inline GPRO::CellSpace<elemType>::
CellSpace(const SpaceDims& dims,
          const elemType &initVal)
{
  if(!initMem(dims) ||
     !initVals(initVal))
  {
    exit(-1);
  }
}

template <class elemType>
inline GPRO::CellSpace<elemType>::
CellSpace(int nRows, int nCols) 
{  
  if(!initMem(SpaceDims(nRows, nCols))) 
  {
    exit(-1); 
  }
}

template <class elemType>
inline GPRO::CellSpace<elemType>::
CellSpace(int nRows, int nCols,
          const elemType &initVal)
{
  if(!initMem(SpaceDims(nRows,nCols)) ||
     !initVals(initVal)) 
  {
    exit(-1);    
  }
}

template <class elemType>
inline GPRO::CellSpace<elemType>::
CellSpace(const CellSpace<elemType> &rhs) 
{
  if(!initMem(rhs._dims)) 
  {
    exit(-1);
  }
  memcpy((void *)_matrix, (const void *)rhs._matrix, 
         sizeof(elemType)*size());
}

template <class elemType>
inline GPRO::CellSpace<elemType>::
~CellSpace() 
{
  clear();
}

template <class elemType> template <class elemType2>
inline bool GPRO::CellSpace<elemType>::
equalDim(const CellSpace<elemType2> &rhs) const 
{
  return (_dims == rhs._dims);
}

template <class elemType>
inline bool GPRO::CellSpace<elemType>::
empty() const 
{
  return (!_matrix &&
          _dims.isNone());
}

template <class elemType>
inline int GPRO::CellSpace<elemType>::
nRows() const
{
  return _dims.nRows();	
}

template <class elemType>
inline int GPRO::CellSpace<elemType>::
nCols() const
{
  return _dims.nCols();	
}

template <class elemType>
inline const GPRO::SpaceDims& GPRO::CellSpace<elemType>::
dims() const 
{
  return _dims;
}

template <class elemType>
inline int GPRO::CellSpace<elemType>::
size() const 
{
  return _dims.size();	
}

template <class elemType>
inline bool GPRO::CellSpace<elemType>::
validCoord(const CellCoord& coord,
           bool warning) const 
{
  bool valid = true;
  if(!coord.valid(_dims)) 
  {
    valid = false;
    if(warning) 
	{
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: coordinates [" << coord \
           << "] out of CellSpace boundary [" \
           << _dims << "]" \
          << endl;
    }
  }
  return valid;
}

template <class elemType>
inline bool GPRO::CellSpace<elemType>::
validCoord(int iRow, int iCol,
           bool warning) const 
{
  return validCoord(CellCoord(iRow, iCol), warning);
}

template <class elemType>
inline bool GPRO::CellSpace<elemType>::
validIdx(int idx,
         bool warning) const 
{
  bool valid = true;
  if(!_dims.validIdx(idx)) 
  {
    valid = false;
    if(warning) 
	{
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: index [" << idx \
           << "] out of CellSpace size [" \
           << size() << "]" \
           << endl;
    }
  }
  return valid;
}

template <class elemType> 
inline elemType* GPRO::CellSpace<elemType>::
operator[](int iRow) 
{
  elemType *pRowElems = 0;
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: empty CellSpace" << endl;
  }
  else 
  {
    if(iRow < 0 || iRow >= nRows()) {
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: row index [" \
           << iRow << "] out of CellSpace boundary [" \
           << _dims << "]" << endl;
    }
    else
	{
      pRowElems = _matrix + iRow*nCols();
    }
  }
  return pRowElems;
}

template <class elemType> 
inline const elemType* GPRO::CellSpace<elemType>::
operator[](int iRow) const 
{
  const elemType *pRowElems = 0;
  if(empty())
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: empty CellSpace" << endl;
  }
  else {
    if(iRow < 0 || iRow >= nRows()) 
	{
      cerr << __FILE__ << " " << __FUNCTION__ \
           << " Error: row index [" \
           << iRow << "] out of CellSpace boundary [" \
           << _dims << "]" << endl;
    }
    else 
	{
      pRowElems = _matrix + iRow*nCols();
    }
  }
  return pRowElems;	
}

template <class elemType> 
inline GPRO::CellSpace<elemType>& GPRO::CellSpace<elemType>::
operator=(const CellSpace<elemType> &rhs)
{
  if(this != &rhs) 
  {
    if(_matrix &&
       _dims != rhs._dims)
	{
      delete[] _matrix;
      _matrix = 0;
    }
    if(!_matrix) 
	{
      if(!initMem(rhs._dims)) 
	  {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: unable to copy a CellSpace" \
             << endl;
        return *this;
      }
    }
    memcpy((void *)_matrix, (const void *)rhs._matrix, 
           sizeof(elemType)*size());
  }
  return *this;
}


template <class elemType> 
inline int GPRO::CellSpace<elemType>::
coord2idx(const CellCoord& coord) const 
{
  return coord.toIdx(_dims);
}

template <class elemType> 
inline int GPRO::CellSpace<elemType>::
coord2idx(int iRow, int iCol) const 
{
  return coord2idx(CellCoord(iRow, iCol));
}

template <class elemType> 
inline const GPRO::CellCoord GPRO::CellSpace<elemType>::
idx2coord(int idx) const 
{
  return CellCoord(idx, _dims);
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
values(vector<elemType> &vVals) const 
{
  bool done = true;
  vVals.erase(vVals.begin(), vVals.end());
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to get values from an empty CellSpace" << endl;
    done = false;
  }
  else 
  {
    for(int iElem = 0; iElem < size(); iElem++) 
	{
      elemType val = _matrix[iElem];
      if(std::find(vVals.begin(), vVals.end(), val) == vVals.end())
	  {
        vVals.push_back(val);
      }
    } // End of iElem loop
    stable_sort(vVals.begin(), vVals.end());
  }
  return done;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
find(IntVect &vFoundIdxs,
     const elemType &val) const
{
  bool found = false;
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to find a value from an empty CellSpace" << endl;
  }
  else 
  {
    elemType *pElem = _matrix;
    while(pElem != _matrix + size()) 
	{
      pElem = std::find(pElem, _matrix+size(), val);
      if(pElem != _matrix + size()) 
	  {
        vFoundIdxs.push_back(pElem - _matrix);
        if(!found) 
		{
          found = true;
        }
        pElem++;  
      }
    } // End of while loop
  }
  return found;      
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
find(IntVect &vFoundIdxs,
     const vector<elemType> &vVals) const 
{
  bool found = false;
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to find a value from an empty CellSpace" << endl;
  }
  else
  {
    for(int iElem = 0; iElem < size(); iElem++) 
	{
      const elemType &elemVal = _matrix[iElem];
      if(std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) 
	  {
        vFoundIdxs.push_back(iElem);
        if(!found) 
		{
          found = true;  
        }
      }
    }
  }
  return found;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
find(IntVect &vFoundIdxs,
     const vector<elemType> &vVals,
     const CoordBR &rectangle) const 
{
  bool found = false;
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to find a value from an empty CellSpace" << endl;
  }
  else 
  {
    for(int iRow = rectangle.minIRow(); 
          iRow <= rectangle.maxIRow(); 
          iRow++) 
	{
      for(int iCol = rectangle.minICol(); 
              iCol <= rectangle.maxICol(); 
              iCol++) 
	  {
        int iElem = coord2idx(iRow, iCol);
        elemType elemVal = _matrix[iElem];
        if(std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) 
		{
          vFoundIdxs.push_back(iElem);
          if(!found) 
		  {
            found = true;  
          }
        }
      }
    }
  }
  return found;
}
     
template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType &val) const 
{
  return std::count(_matrix, _matrix+size(), val);
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType &val,
      const CoordBR &rectangle) const
{
  int found = 0;
  
  if(!rectangle.valid(_dims)) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid bounding rectangle" << endl;
    return found;
  }
  
  for(int iRow = rectangle.minIRow(); 
          iRow <= rectangle.maxIRow(); 
          iRow++) 
  {
    int iColStart = rectangle.minICol(); 
    int iColEnd = rectangle.maxICol();
    int iElemStart = coord2idx(iRow, iColStart);
    int iElemEnd = coord2idx(iRow, iColEnd) + 1;
    found += std::count(_matrix+iElemStart,
                        _matrix+iElemEnd,
                        val);
  } // End of iRow loop
  return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType &val,
      const IntVect &vExcldIdxs) const 
{
  int found = 0;
  for(int iElem = 0; iElem < size(); iElem++) 
  {
    if(std::find(vExcldIdxs.begin(),
                 vExcldIdxs.end(),
                 iElem) == vExcldIdxs.end()) 
	{
      if(_matrix[iElem] == val) 
	  {
        found++;
      }
    }
  } // End of iElem loop
  return found;      
}
                
template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType &val,
      const CoordBR &rectangle,
      const IntVect &vExcldIdxs) const 
{
  int found = 0;
  
  if(!rectangle.valid(_dims)) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid bounding rectangle" << endl;
    return found;
  }
  
  for(int iRow = rectangle.minIRow(); 
          iRow <= rectangle.maxIRow(); 
          iRow++) 
  {
    for(int iCol = rectangle.minICol(); 
            iCol <= rectangle.maxICol(); 
            iCol++) 
	{
      int iElem = coord2idx(iRow, iCol);
      if(std::find(vExcldIdxs.begin(),
                   vExcldIdxs.end(),
                   iElem) == vExcldIdxs.end()) 
	  {
        if(_matrix[iElem] == val) 
		{
          found++;
        }
      }
    } // End of iCol loop
  } // End of iRow loop
  return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType> &vVals) const 
{
  int found = 0;
  for(int iElem = 0; iElem < size(); iElem++) 
  {
    elemType elemVal = _matrix[iElem];
    if(std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) 
	{
      found++;
    }
  }
  return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType> &vVals,
      const CoordBR &rectangle) const 
{
  int found = 0;
  
  if(!rectangle.valid(_dims))
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid bounding rectangle" << endl;
    return found;
  }
  
  for(int iRow = rectangle.minIRow(); 
          iRow <= rectangle.maxIRow(); 
          iRow++)
  {
    for(int iCol = rectangle.minICol(); 
            iCol <= rectangle.maxICol(); 
            iCol++) 
	{
      int iElem = coord2idx(iRow, iCol);
      elemType elemVal = _matrix[iElem];
      if(std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) 
	  {
        found++;
      }
    }
  } // End of iRow loop
  return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType> &vVals,
      const IntVect &vExcldIdxs) const 
{
  int found = 0;
  for(int iElem = 0; iElem < size(); iElem++) 
  {
    if(std::find(vExcldIdxs.begin(),
                 vExcldIdxs.end(),
                 iElem) == vExcldIdxs.end())
	{
      elemType elemVal = _matrix[iElem];
      if(std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) 
	  {
        found++;
      }
    }
  } // End of iElem loop
  return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType> &vVals,
      const CoordBR &rectangle,
      const IntVect &vExcldIdxs) const 
{
  int found = 0;
  
  if(!rectangle.valid(_dims)) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid bounding rectangle" << endl;
    return found;
  }
  
  for(int iRow = rectangle.minIRow(); 
          iRow <= rectangle.maxIRow(); 
          iRow++) 
  {
    for(int iCol = rectangle.minICol(); 
            iCol <= rectangle.maxICol(); 
            iCol++) 
	{
      int iElem = coord2idx(iRow, iCol);
      if(std::find(vExcldIdxs.begin(),
                   vExcldIdxs.end(),
                   iElem) == vExcldIdxs.end()) 
	  {
        elemType elemVal = _matrix[iElem];
        if(std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) 
		{
          found++;
        }
      }
    } // End of iCol loop
  } // End of iRow loop
  return found;
}

template <class elemType> template <class Predicate>
int GPRO::CellSpace<elemType>::
count(Predicate pred) const 
{
  return std::count_if(_matrix, _matrix + size(), pred);
}

template <class elemType> template <class Predicate>
int GPRO::CellSpace<elemType>::
count(Predicate pred,
      const CoordBR &rectangle) const 
{
  int num = 0;
  int iStart, iEnd;
  for(int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) 
  {
    iStart = coord2idx(iRow, rectangle.minICol());
    iEnd = coord2idx(iRow, rectangle.maxICol());
    num += std::count_if(_matrix + iStart,
                         _matrix + iEnd + 1,
                         pred);
  }
  return num;
}

/************************
* iDim = 1 update a row *
* iDim = 2 update a col *
*************************/
template <class elemType>
bool GPRO::CellSpace<elemType>::
updateRowCol(const vector<elemType> &vNewVals,
             int iDim,
             int iRowCol) 
{
  if(empty())
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update an empty CellSpace" << endl;
    return false;
  }

  switch(iDim) 
  {
    case 1: // Update a Row
      if(iRowCol < 0 || iRowCol >= nRows())
	  {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: invalid row index to update a CellSpace" << endl;
        return false;
      }
      if(vNewVals.size() != nCols()) 
	  {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: the size of the new value vector does NOT"
             << " match the number of columns of the CellSpace" << endl;
        return false;
      }
      memcpy((void *)(_matrix+iRowCol*nCols()),
             (const void *)(&vNewVals[0]),
             sizeof(elemType)*nCols());
      break;
    case 2: // Update a Col
      if(iRowCol < 0 || iRowCol >= nCols())
	  {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << "Error: invalid col index to update a CellSpace" << endl;
        return false;
      }
      if(vNewVals.size() != nRows())
	  {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: the size of the new value vector does NOT"
             << " match the number of rows of the CellSpace" << endl;
        return false;
      }
      for(int iRow = 0; iRow < nRows(); iRow++)
	  {
        int iElem = coord2idx(iRow, iRowCol);
        _matrix[iElem] = vNewVals[iRow];
      }
      break;
    default:
      cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalide iDim value" << endl;
      return false;
  }

  return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
updateCells(const vector< GPRO::BasicCell<elemType> > &vCells)
{
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update an empty CellSpace" << endl;
    return false;
  }

  for(int iCell = 0; iCell < vCells.size(); iCell++) 
  {
    CellCoord &coord = vCells[iCell].coord();
    if(!validCoord(coord)) 
	{
      return false;
    }
    int iElem = coord2idx(coord);
    _matrix[iElem] = vCells[iCell].val();
  }

  return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
updateCells(const vector< pair<int, elemType> > &vCells) 
{
  if(empty())
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update an empty CellSpace" << endl;
    return false;
  }

  for(int iCell = 0; iCell < vCells.size(); iCell++) 
  {
    int iElem = vCells[iCell].first;
    if(!validIdx(iElem)) 
	{
      return false;
    }
    _matrix[iElem] = vCells[iCell].second;
  }

  return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
update(Transition<elemType> *pTransition,
       Neighborhood<elemType> *pNbrhood,
       const CoordBR *const pWorkBR) 
{
  if(!pTransition) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update a CellSpace" \
         << " using a void pointer to Transition" \
         << endl;
    return false;
  }
  if(empty()) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update an empty CellSpace" << endl;
    return false;
  }

  const CoordBR *pWBR = pWorkBR;
  CoordBR workBR;
  if(!pWorkBR) 
  {
    if(!pNbrhood->calcWorkBR(workBR, _dims)) 
	{
      return false;
    }
    pWBR = &workBR;
  }
  if(!pWBR->valid(dims()))
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid workBR (" \
         << *pWBR << ")" << endl;
    return false;
  }
  
  if(!pTransition->cellSpace(this) ||
     !pTransition->nbrhood(pNbrhood)) 
  {
    return false;
  }
  
  vector<pair<int, elemType> > vUpdtedCells;
  for(int iRow = pWBR->minIRow();
          iRow <= pWBR->maxIRow(); iRow++) 
  {
    for(int iCol = pWBR->minICol();
            iCol <= pWBR->maxICol(); iCol++) 
	{
      vUpdtedCells.clear();
      CellCoord coord(iRow, iCol);
      if(!pTransition->evaluate(vUpdtedCells, coord)) 
	  {
        return false;
      }
      if(!vUpdtedCells.empty())
	  {
        for(int iCell = 0; iCell < vUpdtedCells.size(); iCell++)
		{
          int &iElem = vUpdtedCells[iCell].first;
          if(!validIdx(iElem)) 
		  {
            cerr << __FILE__ << " " << __FUNCTION__ \
                 << " Error: unable to update coord[" \
                 << coord << "]" << endl;
            return false;
          }
          elemType &val = vUpdtedCells[iCell].second;
          if(pTransition->needFinalize()) 
		  {
            _mUpdtCells[val].push_back(iElem);
          }
          else
		  {
            _matrix[iElem] = val;
          }
        }
      }
    }
  }
  return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
update(IntVect &vIdxs2Eval,
       Transition<elemType> *pTransition,
       Neighborhood<elemType> *pNbrhood,
       const CoordBR *const pWorkBR)
{
  if(!pTransition)
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update a CellSpace" \
         << " using a void pointer to Transition" \
         << endl;
    return false;
  }
  if(empty()) {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to update an empty CellSpace" << endl;
    return false;
  }

  const CoordBR *pWBR = pWorkBR;
  CoordBR workBR;
  if(!pWorkBR)
  {
    if(!pNbrhood->calcWorkBR(workBR, _dims)) 
	{
      return false;
    }
    pWBR = &workBR;
  }
  if(!pWBR->valid(dims()))
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: invalid workBR (" \
         << *pWBR << ")" << endl;
    return false;
  }
  if(!pTransition->cellSpace(this) ||
     !pTransition->nbrhood(pNbrhood))
  {
    return false;
  }

  if(vIdxs2Eval.empty()) 
  {
    return true;
  }
  vector<pair<int, elemType> > vUpdtedCells;
  IntVctItr iIdx = vIdxs2Eval.begin();
  while(iIdx != vIdxs2Eval.end()) 
  {
    int iElem = *iIdx;
    CellCoord coord(iElem, _dims);
    vUpdtedCells.clear();
    if(pWBR->contain(coord))
	{
      if(!pTransition->evaluate(vUpdtedCells, coord))
	  {
        return false;
      }
      if(!vUpdtedCells.empty()) 
	  {
        for(int iCell = 0; iCell < vUpdtedCells.size(); iCell++)
		{
          int &iElem = vUpdtedCells[iCell].first;
          if(!validIdx(iElem)) 
		  {
            cerr << __FILE__ << " " << __FUNCTION__ \
                 << " Error: unable to update coord[" \
                 << coord << "]" << endl;
            return false;
          }
          elemType &val = vUpdtedCells[iCell].second;
          if(pTransition->needFinalize()) 
		  {
            _mUpdtCells[val].push_back(iElem);
          }
          else 
		  {
            _matrix[iElem] = val;
          }
        }
      }
      vIdxs2Eval.erase(iIdx);
    }
    else 
	{
      iIdx++;
    }
  } // End of iIdx loop

  return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
updateFinalize(Transition<elemType> *pTransition,
               Neighborhood<elemType> *pNbrhood)
{
  if(_mUpdtCells.empty()) 
  {
    return true;
  }

  if(pTransition)
  {
    if(!pTransition->cellSpace(this) ||
       !pTransition->nbrhood(pNbrhood))
	{
      return false;
    }
  }

  typename map<elemType, IntVect>::iterator iCellGroup = _mUpdtCells.begin();
  while(iCellGroup != _mUpdtCells.end())
  {
    elemType newVal = (*iCellGroup).first;
    IntVect &vIdxs = (*iCellGroup).second;
    for(int iIdx = 0; iIdx < vIdxs.size(); iIdx++)
	{
      int &idx = vIdxs[iIdx];
      if(!validIdx(idx))
	  {
        cerr << __FILE__ << " " << __FUNCTION__ \
             << " Error: unable to finalize the update of idx[" \
             << idx << "]" << endl;
        return false;
      }
      if(pTransition) 
	  {
        CellCoord coord(idx, _dims);
        if(!pTransition->finalize(newVal, coord)) 
		{
          return false;
        }
      }
      else 
	  {
        _matrix[idx] = newVal;
      }
    }
    iCellGroup++;
  }

  _mUpdtCells.clear();
  return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
add2UpdtMap(const elemType &newVal,
            const CellCoord &coord)
{
  if(!validCoord(coord)) 
  {
    cerr << __FILE__ << " " << __FUNCTION__ \
         << " Error: unable to add a new location to" \
         << " the update cell map" << endl;
    return false;
  }
  _mUpdtCells[newVal].push_back(coord2idx(coord));
  return true;
}

template <class elemType>
void GPRO::CellSpace<elemType>::
cleanUpdtMap() 
{
  _mUpdtCells.clear();
}

#endif
