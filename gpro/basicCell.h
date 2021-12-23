/**
 * \file basicCell
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::BasicCell
 * \version 1.0
 * 
 * \copyright Copyright (c) 2013-2020
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2020 - Wang Yujing - Code reformat
 */

#ifndef BASICCELL_H
#define BASICCELL_H

#include "basicTypes.h"

namespace GPRO {
    /**
     * \ingroup gpro
     * \class BasicCell
     * \brief to record the variable of elemType at the specific coordinate
     */
    template <class elemType>
    class BasicCell {
    public:
        ///< Constructor
        BasicCell() {
        }

        ///< Construct with row&col value
        BasicCell(int iRow, int iCol)
            : _coord(iRow, iCol) {
        }

        ///< Constuct with row&col value and element
        BasicCell(int iRow, int iCol,
                  const elemType& val)
            : _coord(iRow, iCol), _val(val) {
        }

        ///< Construct with coordinate
        BasicCell(const CellCoord& coord)
            : _coord(coord) {
        }

        ///< Construct with coordinate and element
        BasicCell(const CellCoord& coord,
                  const elemType& val)
            : _coord(coord), _val(val) {
        }

        ///< Deep copy
        BasicCell(const BasicCell<elemType>& rhs)
            : _coord(rhs._coord), _val(rhs._val) {
        }

        ///< Destructor
        ~BasicCell() {
        }

        ///< return the copy of the stored value
        const elemType& val() const {
            return _val;
        }

        ///< Getter
        int iRow() const {
            return _coord.iRow();
        }

        ///< Getter
        int iCol() const {
            return _coord.iCol();
        }

        ///< return the copy of coordinate
        const CellCoord& coord() const {
            return _coord;
        }

        ///< Setter
        void val(const elemType& val) {
            _val = val;
        }

        ///< Setter
        void iRow(int iRow) {
            _coord.iRow(iRow);
        }

        ///< Setter
        void iCol(int iCol) {
            _coord.iCol(iCol);
        }

        ///< Shallow copy
        void coord(const CellCoord& coord) {
            _coord = coord;
        }

        ///< Deep copy
        BasicCell& operator=(const BasicCell<elemType>& rhs) {
            if (this != &rhs) {
                _coord = rhs._coord;
                _val = rhs._val;
            }
            return *this;
        }

        bool operator<(const BasicCell<elemType>& rhs) const {
            return _val < rhs._val;
        }

        bool operator<=(const BasicCell<elemType>& rhs) const {
            return _val <= rhs._val;
        }

        bool operator>(const BasicCell<elemType>& rhs) const {
            return _val > rhs._val;
        }

        bool operator>=(const BasicCell<elemType>& rhs) const {
            return _val >= rhs._val;
        }

        bool operator==(const BasicCell<elemType>& rhs) const {
            return _val == rhs._val;
        }

        bool operator!=(const BasicCell<elemType>& rhs) const {
            return _val != rhs._val;
        }

    protected:
        CellCoord _coord; /// coordinate of this cell, including row and col value
        elemType _val; /// value in this cell
    };

};

#endif
