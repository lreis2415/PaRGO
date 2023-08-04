#include"scaOperator.h"

#define Eps 0.000001

void SCAOperator::
d8Layer(RasterLayer<double>& layerD) {
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
void SCAOperator::wslayer(RasterLayer<double>& layerD) {
    _pwsLayer = &layerD;
    Configure(_pwsLayer, false);
}

void SCAOperator::scaLayer(RasterLayer<double>& layerD) {
    _pSCALayer = &layerD;
    Configure(_pSCALayer, true);
}
void SCAOperator::outLayer(RasterLayer<double>& layerD) {
    _pOutLayer = &layerD;
    Configure(_pOutLayer, false);
}
void SCAOperator::initsca(RasterLayer<double>& layerD){
	_pSCALayer=&layerD;
	CellSpace<double>& sca = *(_pSCALayer->cellSpace());
	
	int maxRow = _pSCALayer->_pMetaData->_localworkBR.maxIRow();
    int maxCol = _pSCALayer->_pMetaData->_localworkBR.maxICol();
	int minRow = _pSCALayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pSCALayer->_pMetaData->_localworkBR.minICol();
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			sca[i][j]=0;
		}	
	}
}
void SCAOperator::writesca(RasterLayer<double>& layerD,CellCoord nw){
	_pwkLayer=&layerD;
	CellSpace<double>& wk = *(_pwkLayer->cellSpace());
	CellSpace<double>& sca = *(_pSCALayer->cellSpace());
	//CellSpace<double>& dir = *(_pD8Layer->cellSpace());
	_maxRow = _pwkLayer->_pMetaData->_localworkBR.maxIRow();
    _maxCol = _pwkLayer->_pMetaData->_localworkBR.maxICol();
	int d1[9]={-1,-1,-1,0,0,0,1,1,1};
	int d2[9]={-1,0,1,-1,0,1,-1,0,1};
	int minRow = _pwkLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pwkLayer->_pMetaData->_localworkBR.minICol();
	int gi = nw.iRow()-1;//_pwkLayer->_pMetaData->_MBR.minIRow();
	int gj = nw.iCol();//_pwkLayer->_pMetaData->_MBR.minICol();
	//cout<<"globali:"<<gi<<" globalj:"<<gj<<endl;
	//cout<<"scalayer"<<_pSCALayer->_pMetaData->_localworkBR.maxIRow()<<endl;
	int ti,tj;
	short d;
	double s=0.;
	for(int i=minRow;i<=_maxRow;i++){
		for(int j=minCol;j<=_maxCol;j++){
			//cout<<"sca:"<<sca[i+gi][j+gj]<<endl;
			//cout<<"wk:"<<wk[i][j]<<endl;
			if(fabs(sca[i+gi][j+gj]-_noData)<=Eps||sca[i+gi][j+gj]<wk[i][j]){
				sca[i+gi][j+gj]=wk[i][j];
				//cout<<"wk:"<<wk[i][j]<<" sca:"<<sca[i+gi][j+gj]<<endl;
			}
			
				//cout<<"false"<<endl;

		}
	}
}
void SCAOperator::getarea(int &minrow,int &mincol,int &nrow,int &ncol, int _g,int buf,int id) {
	CellSpace<double>& ws = *(_pwsLayer->cellSpace());

	int pi,pj,qi,qj;
	int minRow = _pwsLayer->_pMetaData->_MBR.minIRow();
	int minCol = _pwsLayer->_pMetaData->_MBR.minICol();
	int maxRow = _pwsLayer->_pMetaData->_MBR.maxIRow();
	int maxCol = _pwsLayer->_pMetaData->_MBR.maxICol();

	pi=-2;
	pj=-2;
	qi=-2;
	qj=-2;
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			if(ws[i][j]==id){
				if(i<pi||pi==-2)
					pi=i;
				if(j<pj||pj==-2)
					pj=j;
				if(i>qi||qi==-2)
					qi=i;
				if(j>qj||qj==-2)
					qj=j;
			}
		}
	}
	
	if((pi-buf)>minRow)
		pi=pi-buf;
	else
		pi=minRow+1;
	if((pj-buf)>minCol)
		pj=pj-buf;
	else
		pj=minCol+1;
	if((qi+buf)<maxRow)
		qi=qi+buf;
	else
		qi=maxRow-1;
	if((qj+buf)<maxCol)
		qj=qj+buf;
	else
		qj=maxCol-1;

	minrow=pi*_g;
	mincol=pj*_g;
	nrow=(qi-pi+1)*_g;
	ncol=(qj-pj+1)*_g;
}
bool SCAOperator::isTermination() {
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

bool SCAOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& d8L = *(_pD8Layer->cellSpace());
    CellSpace<double>& scaL = *(_pSCALayer->cellSpace());
	CellSpace<double>& out = *(_pOutLayer->cellSpace());
    CellSpace<double>& degreeL = *(_degreeLayer.cellSpace());
    Neighborhood<double>& nbrhoodD = *(_pD8Nbrhood);
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    if (num == 0) {

        if (fabs(d8L[iRow][iCol] - _noData) < Eps) {
            degreeL[iRow][iCol] = -2; //init
            scaL[iRow][iCol] = _noData;
		//
        }else if(out[iRow][iCol]>0){
			scaL[iRow][iCol]=0;
			degreeL[iRow][iCol] = 0;
		}
		//
        else {
            degreeL[iRow][iCol] = 0; //init
            scaL[iRow][iCol] = 1;
        }

        if (iRow == _maxRow && iCol == _maxCol) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 1;
            Termination = 0;
        }
        return true;
    }
    if (num == 1) {
        int dir = 8;

        for (int tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
            for (int tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
                if (d8L[tRow][tCol] == dir && (dir + d8L[iRow][iCol] != 8) && dir != 4 && fabs(d8L[iRow][iCol] - _noData) > Eps) {
					if(scaL[tRow][tCol]>=0)
						degreeL[iRow][iCol]++;
                }
                dir--;
            }
		}

		if(degreeL[iRow][iCol] == 0&&out[iRow][iCol]>0){
			scaL[iRow][iCol] += out[iRow][iCol];
		}

        if (iRow == _maxRow && iCol == _maxCol) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 2;
            Termination = 0;
        }
        return true;
    }

    if (degreeL[iRow][iCol] <= 0 && !(iRow == _maxRow && iCol == _maxCol)) {
        return true;
    }

    int dir = 8;
	double s=0.;
    for (int tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
        for (int tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
			if(d8L[tRow][tCol] == dir&&out[tRow][tCol]>0&&out[iRow][iCol]>0)
				s+=out[tRow][tCol];
            if ((degreeL[tRow][tCol] == 0 || degreeL[tRow][tCol] == -1) && d8L[tRow][tCol] == dir) {
                if ((dir + d8L[iRow][iCol] != 8) && dir != 4) {

					//scaL[iRow][iCol] += scaL[tRow][tCol];
					if(out[iRow][iCol]==0)
						scaL[iRow][iCol] += scaL[tRow][tCol];
					else
						scaL[iRow][iCol]+=(scaL[tRow][tCol]-out[tRow][tCol]);
					//
                    degreeL[iRow][iCol]--;
                    if (degreeL[iRow][iCol] == 0) {
                        degreeL[iRow][iCol] = -1;
					}
                }
                if (degreeL[tRow][tCol] == -1) {
                    degreeL[tRow][tCol] = -2;
                }
            }
            dir--;
        }
    }
	//
	if (degreeL[iRow][iCol] == -1) {
			if(s>out[iRow][iCol]&&out[iRow][iCol]>0)
				out[iRow][iCol]=s+1;
			scaL[iRow][iCol] += out[iRow][iCol];
     }
	//
    if (iRow == _maxRow && iCol == _maxCol) {
        MPI_Barrier(MPI_COMM_WORLD);
        int minRow = _pD8Layer->_pMetaData->_localworkBR.minIRow();
        int minCol = _pD8Layer->_pMetaData->_localworkBR.minICol();
		//bool flag=true;
        for (int i = minRow; i <= _maxRow; ++i) {
            for (int j = minCol; j <= _maxCol; ++j) {
                if (degreeL[i][j] == 0) {
                    degreeL[i][j] = -2; //-2 means ending cal., -1 means finishing cal. this itermination
					//flag=false;
					Termination = 0;//add Termination = 0;					
                }
                else {
                    if (degreeL[i][j] == -1) {
                        degreeL[i][j] = 0;
						//flag=false;
						Termination = 0;//add Termination = 0;			
                    }
                }
            }			
        }
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
