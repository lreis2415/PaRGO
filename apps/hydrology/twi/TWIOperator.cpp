#include"TWIOperator.h"

#define Eps 0.0000001

void TWIOperator::SCALayer(RasterLayer<double> &layerD) 
{
	_pSCALayer = &layerD;
	Configure(_pSCALayer, false);

	_cellSize = _pSCALayer->_pMetaData->cellSize;
	_noData = _pSCALayer->_pMetaData->noData;
	_pNbrhood = layerD.nbrhood();
	Neighborhood<double>& nbrhoodD = *(_pNbrhood);
	_iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;

}

void TWIOperator::slopeLayer(RasterLayer<double> &layerD) 
{
	_pSlopeLayer = &layerD;
	Configure(_pSlopeLayer,false);
}

void TWIOperator::twiLayer(RasterLayer<double> &layerD) 
{
	_pTWILayer = &layerD;
	Configure(_pTWILayer,false);
}

bool TWIOperator::isTermination()
{
	return flag;
}

bool TWIOperator::Operator(const CellCoord &coord, bool Operflag)
{
	
	CellSpace<double> &scaL = *(_pSCALayer->cellSpace());
	CellSpace<double> &slopeL = *(_pSlopeLayer->cellSpace());
	CellSpace<double> &twiL = *(_pTWILayer->cellSpace());
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	
	bool loopTag = true;
	for(int tRow = iRow - _iNeighborCells; tRow <= iRow + _iNeighborCells; tRow++)
	{
		for(int tCol = iCol - _iNeighborCells; tCol <= iCol + _iNeighborCells; tCol++)
		{
			if( (scaL[tRow][tCol] == _noData) || (slopeL[tRow][tCol] == _noData) )
			{
				twiL[iRow][iCol] = _noData;
				loopTag = false;
				break;
			}
		}
		if( !loopTag ){
			break;
		}
	}
	if( loopTag ){
		if( slopeL[iRow][iCol]<Eps ){
			slopeL[iRow][iCol] = Eps;
		}

		twiL[iRow][iCol]=log(scaL[iRow][iCol]/(slopeL[iRow][iCol]));
	}

	return true;
}