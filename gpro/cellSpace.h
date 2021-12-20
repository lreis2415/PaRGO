/**
 * \file basicTypes
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::Transition, and GPRO::CellSpace
 * \version 1.0
 * 
 * \copyright Copyright (c) 2013-2020
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 *  This file is adapted from https://github.com/HPSCIL/pRPL.
 * 
 * changelog:
 *  - 1. 2020 - Wang Yujing - Code reformat
 */

#ifndef CELLSPACE_H
#define CELLSPACE_H

#include "basicTypes.h"
#include "basicCell.h"
#include "metaData.h"
#include <iostream>
#include "mpi.h"
#include <gdal_priv.h>
#include <iomanip>

using namespace std;

namespace GPRO {
    template <class elemType>
    class CellSpace;

    template <class elemType>
    class Neighborhood;

    /**
     * \ingroup gpro
     * \class Transition
     * \brief 
     */
    template <class elemType>
    class Transition {
    public:
        Transition(bool onlyUpdtCtrCell = true,
                   bool needFinalize = true,
                   bool needExchange = true,
                   bool edgesFirst = true)
            : _pCellSpace(nullptr),
              _pNbrhood(nullptr),
              _onlyUpdtCtrCell(onlyUpdtCtrCell),
              _needFinalize(needFinalize),
              _needExchange(needExchange),
              _edgesFirst(edgesFirst) {
        }

        virtual ~Transition() {
        }

        bool onlyUpdtCtrCell() const {
            return _onlyUpdtCtrCell;
        }

        bool needFinalize() const {
            return _needFinalize;
        }

        bool needExchange() const {
            return _needExchange;
        }

        bool edgesFirst() const {
            return _edgesFirst;
        }

        void needFinalize(bool ndFnlz) {
            _needFinalize = ndFnlz;
        }

        virtual bool cellSpace(CellSpace<elemType>* pCellSpc) {
            _pCellSpace = pCellSpc;
            if (!_pCellSpace) {
                return false;
            }
            return true;
        }

        virtual bool nbrhood(Neighborhood<elemType>* pNbrhood) {
            _pNbrhood = pNbrhood;
            if (!_pNbrhood) {
                return false;
            }
            return true;
        }

    protected:
        CellSpace<elemType>* _pCellSpace;
        Neighborhood<elemType>* _pNbrhood;
        bool _onlyUpdtCtrCell;
        bool _needFinalize;
        bool _needExchange;
        bool _edgesFirst;
    };

    template <class elemType>
    ostream& operator<<(ostream& os, const CellSpace<elemType>& cellSpace);

    template <class elemType>
    istream& operator>>(istream& is, CellSpace<elemType>& cellSpace);

    /**
     * \ingroup gpro
     * \class CellSpace
     * \brief Data (matrix) of cells
     */
    template <class elemType>
    class CellSpace {
    public:
        CellSpace();
        CellSpace(const SpaceDims& dims);
        CellSpace(const SpaceDims& dims, const elemType& initVal);
        CellSpace(int nRows, int nCols);
        CellSpace(int nRows, int nCols, const elemType& initVal);
        CellSpace(const CellSpace<elemType>& rhs);

        ~CellSpace();

        bool initMem(const SpaceDims& dims); /// Initialize the memory of cell data matrix. Only in private use yet.
        bool initVals(const elemType& initVal); /// Initialize the value of cell data matrix within size.
        void clear(); /// delete the cell data matrix

        template <class elemType2>
        bool equalDim(const CellSpace<elemType2>& rhs) const; /// return if the dimension is the same as param[in] rhs
        bool empty() const; /// check if the matrix and dimension is empty
        int nRows() const; /// return the row num of dimension
        int nCols() const; /// return the column num of dimension
        int size() const; /// return the size (area) of cell dimension.
        const SpaceDims& dims() const;
        bool validCoord(const CellCoord& coord, bool warning = true) const;
        bool validCoord(int iRow, int iCol, bool warning = true) const;
        bool validIdx(int idx, bool warning = true) const;

        elemType* operator[](int iRow);
        const elemType* operator[](int iRow) const;

        ///copy memory
        CellSpace<elemType>& operator=(const CellSpace<elemType>& rhs);

        int coord2idx(const CellCoord& coord) const;
        int coord2idx(int iRow, int iCol) const;
        const CellCoord idx2coord(int idx) const;

        bool values(vector<elemType>& vVals) const;

        /**
         * \brief search globally by value of rasters
         * \param[in] val the value to search.
         * \param[out] vFoundIdxs the 1D-indices of rasters whose values equal val
         */
        bool find(IntVect& vFoundIdxs, const elemType& val) const;
        /**
         * \brief search globally by value of rasters
         * \param[in] val the value to search
         * \param[out] vFoundIdxs the 1D-indices of rasters whose values equal val
         */
        bool find(IntVect& vFoundIdxs, const vector<elemType>& vVals) const;
        /**
         * \brief search globally by values of rasters
         * \param[in] vVals the vector of values to search
         * \param[out] vFoundIdxs the 1D-indices of rasters whose values equal val
         */
        bool find(IntVect& vFoundIdxs, const vector<elemType>& vVals, const CoordBR& rectangle) const;
        /**
         * \brief count number of elements with condition
         * \param[in] val the value to search.
         */
        int count(const elemType& val) const;
        /**
         * \brief count number of elements with condition
         * \param[in] val the value to search
         * \param[in] rectangle area to search within
         */
        int count(const elemType& val, const CoordBR& rectangle) const;
        /**
         * \brief count number of elements with condition
         * \param[in] val the value to search
         * \param[in] vExcldIdxs indices to search within
         */
        int count(const elemType& val, const IntVect& vExcldIdxs) const;
        /**
         * \brief count number of elements with condition
         * \param[in] val the value to search
         * \param[in] rectangle area to search within
         * \param[in] vExcldIdxs indices to search within
         */
        int count(const elemType& val, const CoordBR& rectangle, const IntVect& vExcldIdxs) const;

        int count(const vector<elemType>& vVals) const;
        int count(const vector<elemType>& vVals, const CoordBR& rectangle) const;
        int count(const vector<elemType>& vVals, const IntVect& vExcldIdxs) const;
        int count(const vector<elemType>& vVals, const CoordBR& rectangle, const IntVect& vExcldIdxs) const;

        template <class Predicate>
        int count(Predicate pred) const;

        template <class Predicate>
        int count(Predicate pred, const CoordBR& rectangle) const;
        /**
         * \brief update the values of a row or col of data matrix
         * \param[in] val the value to search
         * \param[in] rectangle area to search within
         * \param[in] vExcldIdxs indices to search within
         */
        bool updateRowCol(const vector<elemType>& vNewVals,
                          int iDim,
                          int iRowCol);
        /**
         * \brief update the value of first n cells, where n is the size of vCells
         * \param[in] vCells input values in the form of BasicCell
         */
        bool updateCells(const vector<BasicCell<elemType>>& vCells);
        /**
         * \brief update the value of first n cells, where n is the size of vCells
         * \param[in] vCells input value in the form of std::pair
         */
        bool updateCells(const vector<pair<int, elemType>>& vCells);

    public:
        elemType* _matrix; /// data of cells
    protected:
        SpaceDims _dims;
        map<elemType, IntVect> _mUpdtCells;
    };

};

/****************************************
*             iostream Methods          *
*****************************************/
template <class elemType>
ostream& GPRO::
operator<<(ostream& os, const CellSpace<elemType>& cellSpace) {
    int nRows = cellSpace.nRows();
    int nCols = cellSpace.nCols();
    //os << nRows << " " << nCols << endl;
    for (int iRow = 0; iRow < nRows; iRow++) {
        os << iRow << ": ";
        for (int iCol = 0; iCol < nCols; iCol++) {
            os << cellSpace[iRow][iCol] << " ";
        } // End of iCol loop
        os << endl;
    } // End of iRow loop
    return os;
}


template <class elemType>
istream& GPRO::
operator>>(istream& is, CellSpace<elemType>& cellSpace) {
    int nRows, nCols;
    is >> nRows;
    is >> nCols;
    SpaceDims dims(nRows, nCols);

    if (!cellSpace.empty() && cellSpace.dims() != dims) {
        cellSpace.clear();
        cellSpace.initMem(dims);
    }
    else if (cellSpace.empty()) {
        cellSpace.initMem(dims);
    }

    for (int iRow = 0; iRow < nRows; iRow++) {
        for (int iCol = 0; iCol < nCols; iCol++) {
            //int iElem = iRow*nCols + iCol;
            is >> cellSpace[iRow][iCol];
        } // End of iCol loop
    } // End of iRow loop
    return is;
}


template <class elemType>
bool GPRO::CellSpace<elemType>::
initMem(const SpaceDims& dims) {
    bool done = true;
    if (!dims.valid()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: invalid Dimensions [" << dims
            << "] to allocate memory for the CellSpace"
            << endl;
        done = false;
    }
    else {
        _dims = dims;
        _matrix = new elemType[dims.size()];
        if (!_matrix) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: unable to allocate memory for the CellSpace"
                << endl;
            done = false;
        }
    }
    return done;
}


template <class elemType>
bool GPRO::CellSpace<elemType>::
initVals(const elemType& initVal) {
    bool done = true;
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to initialize an empty CellSpace"
            << endl;
        done = false;
    }
    for (int iElem = 0; iElem < size(); iElem++) {
        _matrix[iElem] = initVal;
    }
    return done;
}


template <class elemType>
void GPRO::CellSpace<elemType>::
clear() {
    if (_matrix) {
        delete[] _matrix;
        _matrix = 0;
    }
}


template <class elemType>
GPRO::CellSpace<elemType>::
CellSpace()
    : _matrix(nullptr) {
}

template <class elemType>
GPRO::CellSpace<elemType>::
CellSpace(const SpaceDims& dims) {
    if (!initMem(dims)) {
        exit(-1);
    }
}

template <class elemType>
GPRO::CellSpace<elemType>::
CellSpace(const SpaceDims& dims, const elemType& initVal) {
    if (!initMem(dims) || !initVals(initVal)) {
        exit(-1);
    }
}

template <class elemType>
GPRO::CellSpace<elemType>::
CellSpace(int nRows, int nCols) {
    if (!initMem(SpaceDims(nRows, nCols))) {
        exit(-1);
    }
}

template <class elemType>
GPRO::CellSpace<elemType>::
CellSpace(int nRows, int nCols, const elemType& initVal) {
    if (!initMem(SpaceDims(nRows, nCols)) || !initVals(initVal)) {
        exit(-1);
    }
}

template <class elemType>
GPRO::CellSpace<elemType>::
CellSpace(const CellSpace<elemType>& rhs) {
    if (!initMem(rhs._dims)) {
        exit(-1);
    }
    memcpy(static_cast<void*>(_matrix), static_cast<const void*>(rhs._matrix), sizeof(elemType) * size());
}


template <class elemType>
GPRO::CellSpace<elemType>::
~CellSpace() {
    clear();
}


template <class elemType>
template <class elemType2>
bool GPRO::CellSpace<elemType>::
equalDim(const CellSpace<elemType2>& rhs) const {
    return (_dims == rhs._dims);
}


template <class elemType>
bool GPRO::CellSpace<elemType>::
empty() const {
    return (!_matrix && _dims.isNone());
}

template <class elemType>
int GPRO::CellSpace<elemType>::
nRows() const {
    return _dims.nRows();
}

template <class elemType>
int GPRO::CellSpace<elemType>::
nCols() const {
    return _dims.nCols();
}

template <class elemType>
const GPRO::SpaceDims& GPRO::CellSpace<elemType>::
dims() const {
    return _dims;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
size() const {
    return _dims.size();
}


template <class elemType>
bool GPRO::CellSpace<elemType>::
validCoord(const CellCoord& coord, bool warning) const {
    bool valid = true;
    if (!coord.valid(_dims)) {
        valid = false;
        if (warning) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: coordinates [" << coord
                << "] out of CellSpace boundary ["
                << _dims << "]"
                << endl;
        }
    }
    return valid;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
validCoord(int iRow, int iCol, bool warning) const {
    return validCoord(CellCoord(iRow, iCol), warning);
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
validIdx(int idx, bool warning) const {
    bool valid = true;
    if (!_dims.validIdx(idx)) {
        valid = false;
        if (warning) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: index [" << idx
                << "] out of CellSpace size ["
                << size() << "]"
                << endl;
        }
    }
    return valid;
}


template <class elemType>
elemType* GPRO::CellSpace<elemType>::
operator[](int iRow) {
    elemType* pRowElems = nullptr;
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: empty CellSpace" << endl;
    }
    else {
        if (iRow < 0 || iRow >= nRows()) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: row index ["
                << iRow << "] out of CellSpace boundary ["
                << _dims << "]" << endl;
        }
        else {
            pRowElems = _matrix + iRow * nCols();
        }
    }
    return pRowElems;
}

template <class elemType>
const elemType* GPRO::CellSpace<elemType>::
operator[](int iRow) const {
    const elemType* pRowElems = nullptr;
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: empty CellSpace" << endl;
    }
    else {
        if (iRow < 0 || iRow >= nRows()) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: row index ["
                << iRow << "] out of CellSpace boundary ["
                << _dims << "]" << endl;
        }
        else {
            pRowElems = _matrix + iRow * nCols();
        }
    }
    return pRowElems;
}


template <class elemType>
GPRO::CellSpace<elemType>& GPRO::CellSpace<elemType>::
operator=(const CellSpace<elemType>& rhs) {
    if (this != &rhs) {
        if (_matrix &&
            _dims != rhs._dims) {
            delete[] _matrix;
            _matrix = 0;
        }
        if (!_matrix) {
            if (!initMem(rhs._dims)) {
                cerr << __FILE__ << " " << __FUNCTION__
                    << " Error: unable to copy a CellSpace"
                    << endl;
                return *this;
            }
        }
        memcpy(static_cast<void*>(_matrix), static_cast<const void*>(rhs._matrix),
               sizeof(elemType) * size());
    }
    return *this;
}


template <class elemType>
int GPRO::CellSpace<elemType>::
coord2idx(const CellCoord& coord) const {
    return coord.toIdx(_dims);
}

template <class elemType>
int GPRO::CellSpace<elemType>::
coord2idx(int iRow, int iCol) const {
    return coord2idx(CellCoord(iRow, iCol));
}

template <class elemType>
const GPRO::CellCoord GPRO::CellSpace<elemType>::
idx2coord(int idx) const {
    return CellCoord(idx, _dims);
}


template <class elemType>
bool GPRO::CellSpace<elemType>::
values(vector<elemType>& vVals) const {
    bool done = true;
    vVals.erase(vVals.begin(), vVals.end());
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to get values from an empty CellSpace" << endl;
        done = false;
    }
    else {
        for (int iElem = 0; iElem < size(); iElem++) {
            elemType val = _matrix[iElem];
            if (std::find(vVals.begin(), vVals.end(), val) == vVals.end()) {
                vVals.push_back(val);
            }
        } // End of iElem loop
        stable_sort(vVals.begin(), vVals.end());
    }
    return done;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
find(IntVect& vFoundIdxs, const elemType& val) const {
    bool found = false;
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to find a value from an empty CellSpace" << endl;
    }
    else {
        elemType* pElem = _matrix;
        while (pElem != _matrix + size()) {
            pElem = std::find(pElem, _matrix + size(), val);
            if (pElem != _matrix + size()) {
                vFoundIdxs.push_back(pElem - _matrix);
                if (!found) {
                    found = true;
                }
                ++pElem;
            }
        } // End of while loop
    }
    return found;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
find(IntVect& vFoundIdxs, const vector<elemType>& vVals) const {
    bool found = false;
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to find a value from an empty CellSpace" << endl;
    }
    else {
        for (int iElem = 0; iElem < size(); iElem++) {
            const elemType& elemVal = _matrix[iElem];
            if (std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) {
                vFoundIdxs.push_back(iElem);
                if (!found) {
                    found = true;
                }
            }
        }
    }
    return found;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
find(IntVect& vFoundIdxs, const vector<elemType>& vVals, const CoordBR& rectangle) const {
    bool found = false;
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to find a value from an empty CellSpace" << endl;
    }
    else {
        for (int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) {
            for (int iCol = rectangle.minICol(); iCol <= rectangle.maxICol(); iCol++) {
                int iElem = coord2idx(iRow, iCol);
                elemType elemVal = _matrix[iElem];
                if (std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) {
                    vFoundIdxs.push_back(iElem);
                    if (!found) {
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
count(const elemType& val) const {
    return std::count(_matrix, _matrix + size(), val);
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType& val, const CoordBR& rectangle) const {
    int found = 0;
    if (!rectangle.valid(_dims)) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: invalid bounding rectangle" << endl;
        return found;
    }

    for (int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) {
        int iColStart = rectangle.minICol();
        int iColEnd = rectangle.maxICol();
        int iElemStart = coord2idx(iRow, iColStart);
        int iElemEnd = coord2idx(iRow, iColEnd) + 1;
        found += std::count(_matrix + iElemStart, _matrix + iElemEnd, val);
    } // End of iRow loop
    return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType& val, const IntVect& vExcldIdxs) const {

    int found = 0;
    for (int iElem = 0; iElem < size(); iElem++) {
        if (std::find(vExcldIdxs.begin(), vExcldIdxs.end(), iElem) == vExcldIdxs.end()) {
            if (_matrix[iElem] == val) {
                found++;
            }
        }
    } // End of iElem loop
    return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const elemType& val, const CoordBR& rectangle, const IntVect& vExcldIdxs) const {
    int found = 0;
    if (!rectangle.valid(_dims)) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: invalid bounding rectangle" << endl;
        return found;
    }

    for (int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) {
        for (int iCol = rectangle.minICol(); iCol <= rectangle.maxICol(); iCol++) {
            int iElem = coord2idx(iRow, iCol);
            if (std::find(vExcldIdxs.begin(), vExcldIdxs.end(), iElem) == vExcldIdxs.end()) {
                if (_matrix[iElem] == val) {
                    found++;
                }
            }
        } // End of iCol loop
    } // End of iRow loop
    return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType>& vVals) const {
    int found = 0;
    for (int iElem = 0; iElem < size(); iElem++) {
        elemType elemVal = _matrix[iElem];
        if (std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) {
            found++;
        }
    }
    return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType>& vVals, const CoordBR& rectangle) const {
    int found = 0;
    if (!rectangle.valid(_dims)) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: invalid bounding rectangle" << endl;
        return found;
    }

    for (int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) {
        for (int iCol = rectangle.minICol(); iCol <= rectangle.maxICol(); iCol++) {
            int iElem = coord2idx(iRow, iCol);
            elemType elemVal = _matrix[iElem];
            if (std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) {
                found++;
            }
        }
    } // End of iRow loop
    return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType>& vVals, const IntVect& vExcldIdxs) const {
    int found = 0;
    for (int iElem = 0; iElem < size(); iElem++) {
        if (std::find(vExcldIdxs.begin(), vExcldIdxs.end(), iElem) == vExcldIdxs.end()) {
            elemType elemVal = _matrix[iElem];
            if (std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) {
                found++;
            }
        }
    } // End of iElem loop
    return found;
}

template <class elemType>
int GPRO::CellSpace<elemType>::
count(const vector<elemType>& vVals, const CoordBR& rectangle, const IntVect& vExcldIdxs) const {
    int found = 0;
    if (!rectangle.valid(_dims)) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: invalid bounding rectangle" << endl;
        return found;
    }

    for (int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) {
        for (int iCol = rectangle.minICol(); iCol <= rectangle.maxICol(); iCol++) {
            int iElem = coord2idx(iRow, iCol);
            if (std::find(vExcldIdxs.begin(), vExcldIdxs.end(), iElem) == vExcldIdxs.end()) {
                elemType elemVal = _matrix[iElem];
                if (std::find(vVals.begin(), vVals.end(), elemVal) != vVals.end()) {
                    found++;
                }
            }
        } // End of iCol loop
    } // End of iRow loop
    return found;
}

template <class elemType>
template <class Predicate>
int GPRO::CellSpace<elemType>::
count(Predicate pred) const {
    return std::count_if(_matrix, _matrix + size(), pred);
}

template <class elemType>
template <class Predicate>
int GPRO::CellSpace<elemType>::
count(Predicate pred, const CoordBR& rectangle) const {
    int num = 0;
    int iStart, iEnd;
    for (int iRow = rectangle.minIRow(); iRow <= rectangle.maxIRow(); iRow++) {
        iStart = coord2idx(iRow, rectangle.minICol());
        iEnd = coord2idx(iRow, rectangle.maxICol());
        num += std::count_if(_matrix + iStart, _matrix + iEnd + 1, pred);
    }
    return num;
}


/************************
* iDim = 1 update a row *
* iDim = 2 update a col *
*************************/
template <class elemType>
bool GPRO::CellSpace<elemType>::
updateRowCol(const vector<elemType>& vNewVals, int iDim, int iRowCol) {
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to update an empty CellSpace" << endl;
        return false;
    }

    switch (iDim) {
    case 1: // Update a Row
        if (iRowCol < 0 || iRowCol >= nRows()) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: invalid row index to update a CellSpace" << endl;
            return false;
        }
        if (vNewVals.size() != nCols()) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: the size of the new value vector does NOT"
                << " match the number of columns of the CellSpace" << endl;
            return false;
        }
        memcpy(static_cast<void*>(_matrix + iRowCol * nCols()), static_cast<const void*>(&vNewVals[0]),
               sizeof(elemType) * nCols());
        break;
    case 2: // Update a Col
        if (iRowCol < 0 || iRowCol >= nCols()) {
            cerr << __FILE__ << " " << __FUNCTION__
                << "Error: invalid col index to update a CellSpace" << endl;
            return false;
        }
        if (vNewVals.size() != nRows()) {
            cerr << __FILE__ << " " << __FUNCTION__
                << " Error: the size of the new value vector does NOT"
                << " match the number of rows of the CellSpace" << endl;
            return false;
        }
        for (int iRow = 0; iRow < nRows(); iRow++) {
            int iElem = coord2idx(iRow, iRowCol);
            _matrix[iElem] = vNewVals[iRow];
        }
        break;
    default:
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: invalide iDim value" << endl;
        return false;
    }

    return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
updateCells(const vector<BasicCell<elemType>>& vCells) {
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to update an empty CellSpace" << endl;
        return false;
    }

    for (int iCell = 0; iCell < vCells.size(); iCell++) {
        CellCoord& coord = vCells[iCell].coord();
        if (!validCoord(coord)) {
            return false;
        }
        int iElem = coord2idx(coord);
        _matrix[iElem] = vCells[iCell].val();
    }

    return true;
}

template <class elemType>
bool GPRO::CellSpace<elemType>::
updateCells(const vector<pair<int, elemType>>& vCells) {
    if (empty()) {
        cerr << __FILE__ << " " << __FUNCTION__
            << " Error: unable to update an empty CellSpace" << endl;
        return false;
    }

    for (int iCell = 0; iCell < vCells.size(); iCell++) {
        int iElem = vCells[iCell].first;
        if (!validIdx(iElem)) {
            return false;
        }
        _matrix[iElem] = vCells[iCell].second;
    }

    return true;
}

#endif
