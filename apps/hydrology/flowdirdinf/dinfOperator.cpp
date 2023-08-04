#include"dinfoperator.h"

void DinfOperator::
demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
    _pDEMNbrhood = layerD.nbrhood();
    cellSize = _pDEMLayer->_pMetaData->cellSize;
    noData = _pDEMLayer->_pMetaData->noData;
    Configure(_pDEMLayer, false);
}

void DinfOperator::slpLayer(RasterLayer<double>& layerD) {
    _pslpLayer = &layerD;
    Configure(_pslpLayer, false);
}
void DinfOperator::dinfLayer(RasterLayer<double>& layerD) {
    _pDinfLayer = &layerD;
    Configure(_pDinfLayer, false);
}

bool DinfOperator::isTermination() {
    num--;
    return num > 0;
}
void   VSLOPE(float E0,float E1, float E2,double D1,double D2,double DD,float *S,float *A)
{
	float S1,S2,AD;
	if(D1!=0)
		S1=(E0-E1)/D1;
	if(D2!=0)
		S2=(E1-E2)/D2;

	if(S2==0 && S1==0) *A=0;
	else
		*A= (float) atan2(S2,S1);
	AD= (float) atan2(D2,D1);
	if(*A  <   0.)
	{
		*A=0.;
		*S=S1;
	}
	else if(*A > AD)
	{
		*A=AD;
		*S=(E0-E2)/DD;
	}
	else
		*S= (float) sqrt(S1*S1+S2*S2);
}

bool DinfOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& dinf = *(_pDinfLayer->cellSpace());
	CellSpace<double>& slp = *(_pslpLayer->cellSpace());
    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;

    int iRow = coord.iRow();
    int iCol = coord.iCol();

	int ID2[]= {0,1,2,2,1,1,2,2,1 }; 
	int ID1[]= {0,2,1,1,2,2,1,1,2};
	int I1[] = {0,0,-1,-1,0,0,1,1,0 };
	int I2[] = {0,-1,-1,-1,-1,1,1,1,1};
	int J1[] = {0,1,0,0,-1,-1,0,0,1};
	int J2[] = {0,1,1,-1,-1,-1,-1,1,1};
	float  ANGC[]={0,0.,1.,1.,2.,2.,3.,3.,4.};
	float  ANGF[]={0,1.,-1.,1.,-1.,1.,-1.,1.,-1.};

	double lenx,leny;
	lenx=abs(iNeighborCells)*cellSize;
	leny=abs(iNeighborCells)*cellSize;
	double XYsize[3] = {0,lenx,leny};
	double diag = sqrt(lenx*lenx+leny*leny);
	int con;

	if(fabs(dem[iRow][iCol]-noData)<Eps)
	{
		dinf[iRow][iCol]=noData;
		return true;
	}
	//Check if cell is "contaminated" (neighbors have no data)
	con = 0;
	for(int tRow = iRow - iNeighborCells; tRow <= iRow + iNeighborCells; tRow++) {
		for (int tCol = iCol - iNeighborCells; tCol <= iCol + iNeighborCells; tCol++){
			if(con!=-1){
				if( fabs(dem[tRow][tCol]-noData)<Eps ) con=-1;
			}
		}
	}
					
	if( con == -1 ) dinf[iRow][iCol]=noData;

	else {
		float SK[9];
		float ANGLE[9];
		float SMAX=0.;
		int k;
		int KD=0;
		for(k=1; k<=8; k++)
		{
			VSLOPE(
				dem[iRow][iCol],
				dem[iRow+I1[k]][iCol+J1[k]],
				dem[iRow+I2[k]][iCol+J2[k]],
				XYsize[ID1[k]],
				XYsize[ID2[k]],
				diag,
				&SK[k],
				&ANGLE[k]
			);
		}

		dinf[iRow][iCol]=-1; //USE -1 TO INDICATE DIRECTION NOT YET SET 
		for(k=1; k<=8; k++)
		{
			if(SK[k] >  SMAX)
			{
				SMAX=SK[k];
				KD=k;
			}
		}

		if(KD  > 0)
			dinf[iRow][iCol] = (float) (ANGC[KD]*(PI/2)+ANGF[KD]*ANGLE[KD]) ;
		slp[iRow][iCol]=SMAX;	
	}

    return true;
}
