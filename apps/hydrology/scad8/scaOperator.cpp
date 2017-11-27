#include"scaOperator.h"

#define Eps 0.000001

void SCAOperator::
d8Layer(RasterLayer<double> &layerD) 
{
  _pD8Layer = &layerD;
  _pD8Nbrhood = layerD.nbrhood();
  _cellSize = _pD8Layer->_pMetaData->cellSize;
  _noData = _pD8Layer->_pMetaData->noData;
  _maxRow = _pD8Layer->_pMetaData->_localworkBR.maxIRow();
  _maxCol = _pD8Layer->_pMetaData->_localworkBR.maxICol();

  _degreeLayer.copyLayerInfo(*_pD8Layer);
  Configure(_pD8Layer, false);
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
	CellSpace<double> &d8L = *(_pD8Layer->cellSpace());
	CellSpace<double> &scaL = *(_pSCALayer->cellSpace());
	CellSpace<double> &degreeL = *(_degreeLayer.cellSpace());
	Neighborhood<double>& nbrhoodD = *(_pD8Nbrhood);
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	
	if( num ==0 ){

		if( fabs( d8L[iRow][iCol]-_noData )<Eps  ){
			degreeL[iRow][iCol] = -2;	//init
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
				//int objRow = tRow + ((int)d8L[tRow][tCol])/3 -1;
				//int objCol = tCol + ((int)d8L[tRow][tCol])%3 -1;
				if( d8L[tRow][tCol] == dir && (dir + d8L[iRow][iCol] !=8) && dir!=4 && fabs( d8L[iRow][iCol]-_noData )>Eps ){
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

	if( degreeL[iRow][iCol]<=0 && !(iRow == _maxRow && iCol == _maxCol)){
		return true;
	}

	int dir = 8;
	for( int tRow = iRow-1; tRow <= iRow+1; tRow++ ){
		for( int tCol = iCol-1; tCol <= iCol+1; tCol++ ){
			if( (degreeL[tRow][tCol]==0 || degreeL[tRow][tCol]==-1) && d8L[tRow][tCol] == dir ){
				if( (dir + d8L[iRow][iCol] !=8) && dir!=4 ){
					scaL[iRow][iCol] += scaL[tRow][tCol];
					degreeL[iRow][iCol]--;
					if( degreeL[iRow][iCol]==0 ){
						degreeL[iRow][iCol] = -1;
					}
				}
				if(degreeL[tRow][tCol]==-1  ){
					degreeL[tRow][tCol] = -2;
				}

				Termination = 0;
			}
			dir--;
		}
	}

	if( iRow == _maxRow && iCol == _maxCol ){
		MPI_Barrier(MPI_COMM_WORLD);
		int minRow = _pD8Layer->_pMetaData->_localworkBR.minIRow();
		int minCol = _pD8Layer->_pMetaData->_localworkBR.minICol();
		for( int i=minRow; i<=_maxRow; ++i ){
			for( int j=minCol; j<=_maxCol; ++j ){
				if( degreeL[i][j] == 0 ){
					degreeL[i][j] = -2;	//-2 means ending cal., -1 means finishing cal. this itermination
				}else{
					if( degreeL[i][j] == -1 ){
						degreeL[i][j] = 0;
					}
				}
			}
		}
		//num++;
	}

	return true;
}

//another method, only serial
//bool SCAOperator::Operator(const CellCoord &coord,bool operFlag)
//{
//	CellSpace<double> &d8L = *(_pD8Layer->cellSpace());
//	CellSpace<double> &scaL = *(_pSCALayer->cellSpace());
//	CellSpace<double> &degreeL = *(_degreeLayer.cellSpace());
//	Neighborhood<double>& nbrhoodD = *(_pD8Nbrhood);
//	int r = _pD8Layer->_pMetaData->myrank;
//	int iRow = coord.iRow();
//	int iCol = coord.iCol();
//	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
//
//	if( num ==0 ){
//		degreeL[iRow][iCol] = 0;	//init
//		scaL[iRow][iCol] = 1;
//
//		if( iRow == _maxRow && iCol == _maxCol ){
//			MPI_Barrier(MPI_COMM_WORLD);
//			num = 1;
//			Termination = 0;
//		}
//		return true;
//	}
//
//	if( num==1 ){
//		if( d8L[iRow][iCol]>=0 && d8L[iRow][iCol]<9 && d8L[iRow][iCol]!=4 ){
//			int objRow = iRow + ((int)d8L[iRow][iCol])/3 -1;
//			int objCol = iCol + ((int)d8L[iRow][iCol])%3 -1;
//			if( degreeL[objRow][objCol]!=-1 )
//				degreeL[objRow][objCol]++;
//		}else{
//			degreeL[iRow][iCol] = -1;	//done
//		}
//
//		if( iRow == _maxRow && iCol == _maxCol ){
//			MPI_Barrier(MPI_COMM_WORLD);
//			num = 2;
//			Termination = 0;
//		}
//		return true;
//	}
//
//	if( degreeL[iRow][iCol] == 0 ){
//		int objRow = iRow + ((int)d8L[iRow][iCol])/3 -1;
//		int objCol = iCol + ((int)d8L[iRow][iCol])%3 -1;
//		scaL[objRow][objCol] += scaL[iRow][iCol];
//		degreeL[objRow][objCol]--;
//		degreeL[iRow][iCol] = -1;
//
//		//if( fabs(scaL[iRow][iCol] - _noData) >Eps ){
//		//	scaL[iRow][iCol] = scaL[iRow][iCol] / _cellSize;
//		//}
//		//scaL[iRow][iCol] = degreeL[iRow][iCol];
//		Termination = 0;
//	}
//
//	return true;
//}
