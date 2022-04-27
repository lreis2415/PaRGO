#include"PDoperator.h"

void PDOperator::demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
	_pDEMNbrhood = layerD.nbrhood();
    noData = _pDEMLayer->_pMetaData->noData;
	_maxRow = _pDEMLayer->_pMetaData->_localworkBR.maxIRow();
    _maxCol = _pDEMLayer->_pMetaData->_localworkBR.maxICol();

	_smoLayer.copyLayerInfo(*_pDEMLayer);
    Configure(_pDEMLayer, false);
	Configure(&_smoLayer, true);
}

void PDOperator::ucgLayer(RasterLayer<double>& layerD) {
    _pucgLayer = &layerD;
    Configure(_pucgLayer, false);
}

bool PDOperator::isTermination() {
    //num--;
    //return num > 0;
	return true;
}

bool PDOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& ucg = *(_pucgLayer->cellSpace());
	CellSpace<double>& smo_elev = *(_smoLayer.cellSpace());
    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	//coord=处理单元坐标
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    double d=0.;
    int tag = 0;
	if (num == 0) {//smooth the dem
		int dis = 0;
		float wsum=0;
		double elevsum=0;
		float p[3]={0.4,0.1,0.05};
		for(int iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++) {
			for (int iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++) {
				d = dem[iRow1][iCol1];
				dis =abs(iRow1-iRow)+abs(iCol1-iCol);
				if ( dis> 1) {
					elevsum+=p[2]*d;
					wsum+=p[2];
				}
				else if(dis==1){
					elevsum+=p[1]*d;
					wsum+=p[1];
				}
				else{
					elevsum+=p[0]*d;
					wsum+=p[0];
				}
				if (fabs(d - noData) < Eps) {
					tag = 1;
				}
			}
		}
		if(!tag){
			elevsum=elevsum/wsum;
			smo_elev[iRow][iCol]=elevsum;
		}
		else{
			smo_elev[iRow][iCol]=dem[iRow][iCol];
		}
		if (iRow == _maxRow && iCol == _maxCol) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 1;
            Termination = 0;
        }
		return true;
	}

    
	int tCol,tRow;
	//-- Put smoothed elevations back in elevation grid--						
	double emax=smo_elev[iRow][iCol];
	int flag=1;
	int d1[12]={-1,0,-1,0,1,1,1,1,0,0,-1,-1};
	int d2[12]={0,-1,-1,-1,-1,0,0,1,1,1,1,0};
	tag=0;
	d=0;	
	for (int k=0;k<12;k++){
		tCol=iCol+d1[k];
		tRow=iRow+d2[k];
		d=smo_elev[tRow][tCol];
		if(fabs(d - noData) < Eps)
		{
			tag=1;
		}
		else {		
			if(d>emax)
				flag=0;	
			if(k%3==2)
			{
				if(flag==1)
					tag=1;
				flag=1;
			}
		}
	}
	if(tag==1)
		ucg[iRow][iCol]=0;
	else
		ucg[iRow][iCol]=1;
	return true;

}
	

			
		
		

