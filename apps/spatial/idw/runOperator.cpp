#include "runOperator.h"

void RUNOperator::demLayer(RasterLayer<double> &layerD) 
{
	_pDEMLayer = &layerD;
	_pDEMNbrhood = layerD.nbrhood();
	Configure(_pDEMLayer, false);

	_noData = _pDEMLayer->_pMetaData->noData;
	_cellSize = _pDEMLayer->_pMetaData->cellSize;
	_xSize = _pDEMLayer->_pMetaData->_localdims.nRows();
	_ySize = _pDEMLayer->_pMetaData->_localdims.nCols();
	_nRows = _pDEMLayer->_pMetaData->row;
	_nCols = _pDEMLayer->_pMetaData->column;
	_rank = _pDEMLayer->_pMetaData->myrank;
}

void RUNOperator::pitLayer(RasterLayer<double> &layerD){
	_pPitLayer = &layerD;
	Configure(_pPitLayer, false);
}

bool RUNOperator::isTermination()
{
	return flag;
}

bool RUNOperator::Operator(const CellCoord &coord, bool operFlag)
{
	CellSpace<double> &demL = *(_pDEMLayer->cellSpace());
	CellSpace<double> &pitL = *(_pPitLayer->cellSpace());
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;

	int iRow = coord.iRow();
	int iCol = coord.iCol();

	if((iRow==_xSize-2)&&(iCol==_ySize-2)){
		MPI_Barrier(MPI_COMM_WORLD); 
		num++;
		if(num==0)
			Termination=0;
	}
	if( demL[iRow][iCol]==_noData ){
		pitL[iRow][iCol] = _noData;
		return true;
	}

	for(int i = iRow - iNeighborCells; i <= iRow + iNeighborCells; i++){
		for(int j = iCol - iNeighborCells; j <= iCol + iNeighborCells; j++){
			pitL[iRow][iCol] = demL[iRow][iCol];
		}
	}


	//这一堆代码只是之前为了测试效率加速，并无实际意义
	//if( num==0 ){
	//	//if( iRow==1 && iCol==1 && _pDEMLayer->_pMetaData->myrank==0 )
	//	//	cout<<_noData<<endl;
	//	if(fabs(demL[iRow][iCol]-_noData)>Eps){
	//		//pitL[iRow][iCol]=100000.0;
	//		pitL[iRow][iCol]=demL[iRow][iCol];
	//	}else{
	//		pitL[iRow][iCol] = _noData;
	//	}

	//	if((iRow==_xSize-2)&&(iCol==_ySize-2)){
	//		MPI_Barrier(MPI_COMM_WORLD); 
	//		num++;
	//		Termination=0;
	//	}
	//}else{
	//	double gap = 0.0005;
	//	if (fabs(demL[iRow][iCol]-_noData)>Eps){
	//		if( pitL[iRow][iCol] - demL[iRow][iCol] < 0.001 ){
	//		//if( num>11 ){
	//		//	cout<<_pDEMLayer->_pMetaData->myrank<<" "<<iRow<<" "<<iCol<<endl;
	//		//}
	//			pitL[iRow][iCol] += gap;
	//			//if( iRow==1 && iCol==525 && _pDEMLayer->_pMetaData->myrank==7 ){
	//			//	cout<<_noData<<endl;
	//			//	cout<<pitL[iRow][iCol]<<" "<<demL[iRow][iCol]<<" num="<<num<<endl;
	//			//}
	//			Termination = 0;
	//			double d[9];
	//			int k = 0;
	//			int tag=0;
	//			for( int tmp=0; tmp<20; ++tmp ){
	//				//这个循环里的内容都只是为了增加单栅格计算代价，没别的含义
	//				k = 0;
	//				tag=0;
	//				for(int i = iRow - iNeighborCells; i <= iRow + iNeighborCells; i++){
	//					for(int j = iCol - iNeighborCells; j <= iCol + iNeighborCells; j++){
	//						d[k] = demL[i][j];
	//						if(d[k]==_noData)
	//						{
	//							tag=1;
	//						}
	//						k++;
	//					}
	//				}
	//				if(tag==1)
	//				{
	//					double tmp = _noData;
	//					return true;
	//				}else{
	//					double dx = (d[8] + 2*d[5] + d[2] - d[6] -2*d[3] - d[0])/(8.0*_cellSize);
	//					double dy = (d[2] + 2*d[1] + d[0] - d[6] -2*d[7] - d[8])/(8.0*_cellSize);
	//					double tmp = sqrt(dx*dx + dy*dy);
	//				}
	//			}
	//		}
	//		//if( pitL[iRow][iCol] > demL[iRow][iCol] ){
	//		//	for(int i = iRow - iNeighborCells; i <= iRow + iNeighborCells; i++){
	//		//		for(int j = iCol - iNeighborCells; j <= iCol + iNeighborCells; j++){
	//		//			if( demL[iRow][iCol] >= (pitL[i][j] + gap) ){
	//		//				pitL[iRow][iCol] = demL[iRow][iCol];
	//		//				Termination = 0;
	//		//			}else {
	//		//				if( pitL[iRow][iCol] > (pitL[i][j] + gap) ){
	//		//					pitL[iRow][iCol] = pitL[i][j] + gap;
	//		//					Termination = 0;
	//		//				}
	//		//			}
	//		//		}
	//		//	}
	//		//}
	//	}
	//	if((iRow==_xSize-2)&&(iCol==_ySize-2)){
	//		++num;
	//		//if(num%50==0)
	//		//	cout<<"pitIterNum"<<num<<" Termination"<<Termination<<endl;
	//		MPI_Barrier(MPI_COMM_WORLD); 
	//		int tmpTer = 0;
	//		MPI_Allreduce(&Termination,&tmpTer,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	//		if(tmpTer != _pDEMLayer->_pMetaData->processor_number){
	//			Termination = 0;
	//		}
	//	}
	//}

	return true;
}