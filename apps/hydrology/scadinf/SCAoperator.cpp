#include"dinfoperator.h"

void DinfOperator::
dinfLayer(RasterLayer<double>& layerD) {
    _pDinfLayer = &layerD;
    _pDinfNbrhood = layerD.nbrhood();
    cellSize = _pDinfLayer->_pMetaData->cellSize;
    noData = _pDinfLayer->_pMetaData->noData;
	_maxRow = _pDinfLayer->_pMetaData->_localworkBR.maxIRow();
    _maxCol = _pDinfLayer->_pMetaData->_localworkBR.maxICol();

    _degreeLayer.copyLayerInfo(*_pDinfLayer);
    Configure(_pDinfLayer, false);
    Configure(&_degreeLayer, true);
	

}

void DinfOperator::scaLayer(RasterLayer<double>& layerD) {
    _pSCALayer = &layerD;
    Configure(_pSCALayer, true);
}

bool DinfOperator::isTermination() {
    //num--;
    return true;
}

double prop( double a, int k, double dx1 , double dy1) {

	double aref[10] = { -atan2(dy1,dx1), 0., -aref[0],(double)(0.5*PI),PI-aref[2],(double)PI,
                       PI+aref[2],(double)(1.5*PI),2.*PI-aref[2],(double)(2.*PI) };
	double p=0.;
	if(k<=0)k=k+8;  
	if(k == 1 && a > PI)a=(float)(a-2.0*PI);
	if(a > aref[k-1] && a < aref[k+1]){
		if( a > aref[k])
			p=(aref[k+1]-a)/(aref[k+1]-aref[k]);
		else
			p=(a-aref[k-1])/(aref[k]-aref[k-1]);
	}

	if( p < 1e-5) return -1.;
	else return(p);
}


bool DinfOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dinf = *(_pDinfLayer->cellSpace());
    CellSpace<double>& scaL = *(_pSCALayer->cellSpace());
	CellSpace<double>& degreeL = *(_degreeLayer.cellSpace());
	
    Neighborhood<double>& nbrhoodD = *(_pDinfNbrhood);
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	
    int iRow = coord.iRow();
    int iCol = coord.iCol();
	
	
	short d;
	double dir;
	double p=0.;
	double lenx,leny;
	lenx=iNeighborCells*cellSize;
	leny=iNeighborCells*cellSize;
	int I1[] = {0,-1,-1,0,0,1,1,0 };
	int J1[] = {1,0,0,-1,-1,0,0,1};
	int I2[] = {-1,-1,-1,-1,1,1,1,1};
	int J2[] = {1,1,-1,-1,-1,-1,1,1};
	int d1[] = {0,0,-1,-1,-1,0,1,1,1};
	int d2[] = {0,1,1,0,-1,-1,-1,0,1};
	
	if (num == 0) {

        if (fabs(dinf[iRow][iCol] - noData) < Eps) {
            degreeL[iRow][iCol] = -2; //init
            scaL[iRow][iCol] = noData;
        }
        else {
            degreeL[iRow][iCol] = 0; //init
            scaL[iRow][iCol] = 1.;
        }
        if (iRow == _maxRow && iCol == _maxCol) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 1;
            Termination = 0;
        }
        return true;
    }
	
	int tRow,tCol;
    if (num == 1) {
        
		for(int k=1;k<=8;k++){
			tRow=d1[k]+iRow;
			tCol=d2[k]+iCol;
			if(!fabs(dinf[tRow][tCol]-noData)<Eps){
				dir=dinf[tRow][tCol];
				p=prop(dir,(k+4)%8,lenx,leny);

				if(p>0.0){
				
					if(!fabs(dinf[iRow][iCol]-noData)<Eps){
						degreeL[iRow][iCol]++;
					}
				}
			}
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

	
	
	for(int k=1;k<=8;k++){
			tRow=d1[k]+iRow;
			tCol=d2[k]+iCol;
			dir=dinf[tRow][tCol];
			if(fabs(dir-noData)>Eps&&(degreeL[tRow][tCol]==0)){
				p=prop(dir,(k+4)%8,lenx,leny);

				if(p>0.0){
					scaL[iRow][iCol] =scaL[iRow][iCol] + scaL[tRow][tCol]*p;
					degreeL[iRow][iCol]--;
					
					if (degreeL[iRow][iCol] == 0) 
						degreeL[iRow][iCol] = -1;
					
				}
		}
	}

    if (iRow == _maxRow && iCol == _maxCol) {
		
        MPI_Barrier(MPI_COMM_WORLD);
        int minRow = _pDinfLayer->_pMetaData->_localworkBR.minIRow();
        int minCol = _pDinfLayer->_pMetaData->_localworkBR.minICol();
		for (int i = minRow; i <= _maxRow; ++i) {
            for (int j = minCol; j <= _maxCol; ++j) {
				if (degreeL[i][j] == 0) {
                    degreeL[i][j] = -2; //-2 means ending cal., -1 means finishing cal. this itermination
					Termination = 0;
                }
                else {
                    if (degreeL[i][j] == -1) {
                        degreeL[i][j] = 0;
						Termination = 0;
                    }
                }
			}
		}
	}

    return true;
}
