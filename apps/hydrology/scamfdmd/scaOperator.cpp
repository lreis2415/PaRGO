#include"scaOperator.h"

#define Eps 0.000001

void SCAOperator::
mfdLayer(vector<RasterLayer<double> *> &layerV) 
{
	for( int i=0; i<8; ++i ){
		_weightLayerVec.push_back(layerV[i]);
		Configure(_weightLayerVec[i], false);
		weightLs.push_back(_weightLayerVec[i]->cellSpace());
	}

  _pNbrhood = layerV[0]->nbrhood();
  _cellSize = layerV[0]->_pMetaData->cellSize;
  _noData = layerV[0]->_pMetaData->noData;
  _maxRow = layerV[0]->_pMetaData->_localworkBR.maxIRow();
  _maxCol = layerV[0]->_pMetaData->_localworkBR.maxICol();
  Neighborhood<double>& nbrhoodD = *(_pNbrhood);
  _iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;

  _degreeLayer.copyLayerInfo(*layerV[0]);
  Configure(&_degreeLayer, true);
}

void SCAOperator::scaLayer(RasterLayer<double> &layerD) 
{
  _pSCALayer = &layerD;
  Configure(_pSCALayer,true);
}

bool SCAOperator::isTermination()
{
	return true;
	//num--;
	//if(num > 0)
	//{
	//	return true;
	//}
	//else
	//{
	//	return false;
	//}
}

bool SCAOperator::Operator(const CellCoord &coord,bool operFlag)
{
	CellSpace<double> &scaL = *(_pSCALayer->cellSpace());
	CellSpace<double> &degreeL = *(_degreeLayer.cellSpace());
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	
	if( num ==0 ){
		if( fabs( (*weightLs[0])[iRow][iCol]-_noData )<Eps  ){
			degreeL[iRow][iCol] = -1;	//init
			scaL[iRow][iCol] = _noData;
		}else{
			degreeL[iRow][iCol] = 0;	//init
			scaL[iRow][iCol] = 1;
		}

		if( iRow == _maxRow && iCol == _maxCol ){
			MPI_Barrier(MPI_COMM_WORLD);
			num = 1;
			Termination = 0;
		}
		return true;
	}

	if( num==1 ){
		int dir = 8;
		for( int tRow = iRow-1; tRow <= iRow+1; tRow++ ){
			for( int tCol = iCol-1; tCol <= iCol+1; tCol++ ){
				if( dir>4 && (*weightLs[dir-1])[tRow][tCol]>0 ){
					degreeL[iRow][iCol]++;
				}
				if( dir<4 && (*weightLs[dir])[tRow][tCol]>0 ){
					degreeL[iRow][iCol]++;
				}
				dir--;
			}
		}
		if( iRow == _maxRow && iCol == _maxCol ){
			MPI_Barrier(MPI_COMM_WORLD);
			num = 2;
			Termination = 0;
		}
		return true;
	}

	if( degreeL[iRow][iCol]<=0 ){
		return true;
	}
	int dir = 8;
	for( int tRow = iRow-1; tRow <= iRow+1; tRow++ ){
		for( int tCol = iCol-1; tCol <= iCol+1; tCol++ ){

			if( dir>4 && degreeL[tRow][tCol]==0 && (*weightLs[dir-1])[tRow][tCol]>0 ){
				scaL[iRow][iCol] += scaL[tRow][tCol]*(*weightLs[dir-1])[tRow][tCol];
				degreeL[iRow][iCol]--;
				(*weightLs[dir-1])[tRow][tCol] = 0;
				Termination = 0;
			}
			if( dir<4 && degreeL[tRow][tCol]==0 && (*weightLs[dir])[tRow][tCol]>0 ){
				scaL[iRow][iCol] += scaL[tRow][tCol]*(*weightLs[dir])[tRow][tCol];
				degreeL[iRow][iCol]--;
				(*weightLs[dir])[tRow][tCol] = 0;
				Termination = 0;
			}

			dir--;
		}
	}

	return true;
}
