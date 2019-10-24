#ifndef BASICCELL_H
#define BASICCELL_H

/***************************************************************************
* basicCell.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::BasicCell
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

namespace GPRO 
{
  template <class elemType>
  class BasicCell 
  {
    public:
      // Constructor
      BasicCell() {}

      BasicCell(int iRow, int iCol)
        :_coord(iRow, iCol){}
      
      BasicCell(int iRow, int iCol,
                const elemType &val)
        :_coord(iRow, iCol),
         _val(val) {}
      
      BasicCell(const CellCoord &coord)
        :_coord(coord) {}
      
      BasicCell(const CellCoord &coord,
                const elemType &val)
        :_coord(coord),
         _val(val) {}

      BasicCell(const BasicCell<elemType> &rhs)
        :_coord(rhs._coord),
         _val(rhs._val) {}

      // Deconstructor
      ~BasicCell() {}

      const elemType& val() const 
	  {
        return _val;
      }
      int iRow() const 
	  {
        return _coord.iRow();
      }
      int iCol() const 
	  {
        return _coord.iCol();
      }
      const CellCoord& coord() const 
	  {
        return _coord;
      }

      void val(const elemType &val) 
	  {
        _val = val;
      }
      void iRow(int iRow) 
	  {
        _coord.iRow(iRow);
      }
      void iCol(int iCol) 
	  {
        _coord.iCol(iCol);
      }
      void coord(const CellCoord &coord) 
	  {
        _coord = coord;
      }

      BasicCell& operator=(const BasicCell<elemType> &rhs) 
	  {
        if(this != &rhs) 
		{
          _coord = rhs._coord;
          _val = rhs._val;
        }
        return *this;
      }

      bool operator<(const BasicCell<elemType> &rhs) const 
	  {
        return _val < rhs._val;
      }
      bool operator<=(const BasicCell<elemType> &rhs) const 
	  {
        return _val <= rhs._val;
      }
      bool operator>(const BasicCell<elemType> &rhs) const 
	  {
        return _val > rhs._val;
      }
      bool operator>=(const BasicCell<elemType> &rhs) const 
	  {
        return _val >= rhs._val;
      }
      bool operator==(const BasicCell<elemType> &rhs) const 
	  {
        return _val == rhs._val;
      }
      bool operator!=(const BasicCell<elemType> &rhs) const 
	  {
        return _val != rhs._val;
      }
    protected:
      CellCoord _coord;
      elemType _val;
  };

};

#endif
