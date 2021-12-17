/**
 * \file basicCell
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::WeightedCell
 * \version 1.0
 * 
 * \copyright Copyright (c) 2013
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2020 - Yujing Wang - Code reformat
 */

#ifndef WEIGHTEDCELL_H
#define WEIGHTEDCELL_H

#include "basicTypes.h"
#include "basicCell.h"

namespace GPRO {
    /**
     * \ingroup gpro
     * \class WeightedCell
     * \brief BasicCell with weight, used as neighbor cells
     */
    template <class elemType>
    class WeightedCell : public BasicCell<elemType> {
    public:
        WeightedCell()
            : BasicCell<elemType>(), _weight(0.0) {
        }

        WeightedCell(int iRow, int iCol)
            : BasicCell<elemType>(iRow, iCol), _weight(1.0) {
        }

        WeightedCell(const CellCoord& coord)
            : BasicCell<elemType>(coord), _weight(1.0) {
        }

        WeightedCell(int iRow, int iCol, double weight)
            : BasicCell<elemType>(iRow, iCol), _weight(weight) {
        }

        WeightedCell(const CellCoord& coord, double weight)
            : BasicCell<elemType>(coord), _weight(weight) {
        }

        WeightedCell(int iRow, int iCol, const elemType& val, double weight)
            : BasicCell<elemType>(iRow, iCol, val), _weight(weight) {
        }

        WeightedCell(const CellCoord& coord, const elemType& val, double weight)
            : BasicCell<elemType>(coord, val), _weight(weight) {
        }

        WeightedCell(const WeightedCell<elemType>& rhs)
            : BasicCell<elemType>(BasicCell<elemType>(rhs)), _weight(rhs._weight) {
        }

        WeightedCell(const BasicCell<elemType>& rhs)
            : BasicCell<elemType>(rhs), _weight(1.0) {
        }

        ~WeightedCell() {
        }

        double weight() const {
            return _weight;
        }

        void weight(double weight) {
            _weight = weight;
        }

        WeightedCell& operator=(const WeightedCell<elemType>& rhs) {
            if (this != &rhs) {
                BasicCell<elemType>::operator=(BasicCell<elemType>(rhs));
                _weight = rhs._weight;
            }
            return *this;
        }

    protected:
        double _weight; /// weight of the basic cell
    };

};

#endif
