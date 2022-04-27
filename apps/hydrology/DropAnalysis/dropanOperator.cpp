#include"dropanOperator.h"

void DropanOperator::
ssaLayer(RasterLayer<double>& layerD) {
    _pssaLayer = &layerD;
	_pssaNbrhood = layerD.nbrhood();
	_noData = _pssaLayer->_pMetaData->noData;
	_maxrow = _pssaLayer->_pMetaData->_localworkBR.maxIRow();
    _maxcol = _pssaLayer->_pMetaData->_localworkBR.maxICol();
	_rank=_pssaLayer->_pMetaData->myrank;
	_elevOut.copyLayerInfo(*_pssaLayer);
	_contribs.copyLayerInfo(*_pssaLayer);
	_orderout.copyLayerInfo(*_pssaLayer);

    Configure(_pssaLayer, false);
	Configure(&_elevOut, true);
	Configure(&_contribs, true);
	Configure(&_orderout, true);
}


void DropanOperator::demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
	_pssaNbrhood = layerD.nbrhood();
    Configure(_pDEMLayer, false);
}



void DropanOperator::dirLayer(RasterLayer<double>& layerD) {
    _pDirLayer = &layerD;
	_pssaNbrhood = layerD.nbrhood();
    Configure(_pDirLayer, false);
}
void DropanOperator::areaLayer(RasterLayer<double>& layerD) {
    _pareaLayer = &layerD;
	_pssaNbrhood = layerD.nbrhood();
	_cellSize = _pareaLayer->_pMetaData->cellSize;
    Configure(_pareaLayer, false);
}

bool DropanOperator::isTermination() {
    //num--;
    //return num > 0;
	return true;
}


//does the appropriate updates when a junction is found
void updateAtJunction(short oOut,long i, long ni,long j, long nj, CellSpace<double>& dir,
                      CellSpace<double>& orderout, CellSpace<double>& elevOut, 
                      CellSpace<double>& dem, bool &newstream,
                      double &s1,double &s1sq,double &s2,double &s2sq,long &n1,long &n2){
	short o;
	double drop;

	//get the order of the pointing cell
	o=orderout[ni][nj];
	
	// if the order increases here, there is an end to the stream segment
	if(o<oOut){
		//  This is the terminus of a stream so accumulate drops
		drop = elevOut[ni][nj]-dem[i][j];
		if(o==1){
			s1= s1+drop;
			s1sq= s1sq+drop*drop;
			n1= n1+1;
		}else{
			s2= s2+drop;
			s2sq= s2sq+drop*drop;
			n2= n2+1;
		}
	}else{
		//  This is the continuation of a main stream to pass elevOut on down
		elevOut[i][j]=elevOut[ni][nj];
		newstream=false;
	}

}

//find the orderOut of a cell based on all the incoming orders
short newOrder(short nOrder[9], bool &junction, bool &source){
	short i,j,temp;
	short oOut=1;  //  initialize to 1
	short ordermax=0;
	short count=0;
	junction=false;
	source=true;
	for(i=0; i<9; i++)
	{
		if(nOrder[i]>0)  // Here an inflowing stream
		{
			count=count+1;
			source=false;
			if(count == 1) //  First inflowing stream
			{
				oOut=nOrder[i];
				ordermax=nOrder[i];
			}else  // Here count is > 1
			{
				if(nOrder[i]> oOut)
				{
					ordermax = nOrder[i];
					oOut=nOrder[i];
				}
				else if(nOrder[i]==oOut) oOut=ordermax+1;
				//  This logic to ensure that order is max of highest in or second highest in +1
				junction=true;
			}
		}
	}
	return oOut;
}


bool DropanOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& ssa = *(_pssaLayer->cellSpace());
	CellSpace<double>& elevout = *(_elevOut.cellSpace());
	CellSpace<double>& contribs = *(_contribs.cellSpace());
	CellSpace<double>& orderout = *(_orderout.cellSpace());
	CellSpace<double>& dir = *(_pDirLayer->cellSpace());
	CellSpace<double>& area = *(_pareaLayer->cellSpace());
	Neighborhood<double>& nbrhoodD = *(_pssaNbrhood);

    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	int iRow = coord.iRow();
    int iCol = coord.iCol();
	int d1[9]={-1,-1,-1,0,0,0,1,1,1};
	int d2[9]={-1,0,1,-1,0,1,-1,0,1};
	
	short d;
	int nexti,nextj;
	int minRow = _pDirLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pDirLayer->_pMetaData->_localworkBR.minICol();
	
	
	if(num==0)
	{
		if(fabs(ssa[iRow][iCol] - _noData) > Eps&&ssa[iRow][iCol]>=thresh)
			contribs[iRow][iCol]=0;
		else
			contribs[iRow][iCol]=-2;//-2 means ending cal.,
		
		if (area[iRow][iCol] > maxval) {
			d=dir[iRow][iCol];
			nexti=iRow+d1[d];
			nextj=iCol+d2[d];
			if(fabs(ssa[nexti][nextj]-_noData)<Eps){
				maxval=area[iRow][iCol];
			}
			
			//cout<<"maxval"<<maxval<<endl;
		}
		if (iRow == _maxrow && iCol == _maxcol) {
			MPI_Barrier(MPI_COMM_WORLD);
			double maxv;
			maxv=maxval;
			MPI_Allreduce(&maxv,&totalarea,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			num = 1;
			Termination = 0;
        }
		return true;
	}
	
	int tRow,tCol;
	
	if(num==1){
		
		if(fabs(ssa[iRow][iCol] - _noData) > Eps&&ssa[iRow][iCol]>=thresh){	
			for (tRow = iRow - 1; tRow <= iRow + 1; tRow++) {
				for (tCol = iCol - 1; tCol <= iCol + 1; tCol++) {
					if(!fabs(dir[tRow][tCol]-_noData)<Eps&&dir[tRow][tCol]>=0){
						d=dir[tRow][tCol];
						if (d1[d]+tRow==iRow && d2[d]+tCol==iCol&&d + dir[iRow][iCol] != 8&& d != 4 ) {
						
							if(ssa[tRow][tCol]>=thresh&&fabs(dem[iRow][iCol] - _noData) > Eps)contribs[iRow][iCol]++;

						}
					}

               
				}
			}
		}
			
		if (iRow == _maxrow && iCol == _maxcol) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 2;
            Termination = 0;	
        }
		return true;
	}
	
	if (contribs[iRow][iCol] < 0 && !(iRow == _maxrow && iCol == _maxcol)) {
        return true;
    }


	if(num>=2){
		if(!contribs[iRow][iCol]){
			long i,j,pi,pj,pd,k;
			double lenx,leny;

			short nOrder[9];  // neighborOrders
			//put the order of all the neighboring cells into an array, and save the i,j, and elev from contributor.
			//pi,pj, and pd are useful only if there is no junction, where there is necessarily one and only one contributor.
			for(k=0;k<9;++k)nOrder[k]=0;
			k=0;
			for(tRow = iRow - iNeighborCells; tRow <= iRow + iNeighborCells; tRow++) {
				for (tCol = iCol - iNeighborCells; tCol <= iCol + iNeighborCells; tCol++){
					if(dir[tRow][tCol]>=0){
						d=dir[tRow][tCol];

						if(d1[d]+tRow ==iRow&& d2[d]+tCol==iCol &&d + dir[iRow][iCol] != 8&& d != 4  
							&& !fabs(orderout[tRow][tCol]-_noData)<Eps){				
							nOrder[k]=orderout[tRow][tCol];
							
						//  Accumulate length
							pd=k; 
							pi=tRow;
							pj=tCol;
							lenx=abs(tRow-iRow);
							leny=abs(tCol-iCol);
							length=length+sqrt(lenx*lenx+leny*leny);
						}
					}
					k++;
				}
			}			
		// Determine the order of this cell
			short oOut;
			bool source;
			bool junction; // junction set to true/false in newOrder
			oOut = newOrder(nOrder,junction,source);
			
			orderout[iRow][iCol]=oOut;
			if(source){
				elevout[iRow][iCol]=dem[iRow][iCol];
			}else if(!junction){
			// if not a junction, transfer elevOut
				elevout[iRow][iCol]=elevout[pi][pj];
			}else{

			// if it is a junction, update global values
				bool newstream=true;  // Flag to indicate whether the junction results in a new Strahler stream
				for(tRow = iRow - iNeighborCells; tRow <= iRow + iNeighborCells; tRow++) {
					for (tCol = iCol - iNeighborCells; tCol <= iCol + iNeighborCells; tCol++){
						if(dir[tRow][tCol]>=0){
							d=dir[tRow][tCol];
							if((d1[d]+tRow ==iRow&& d2[d]+tCol==iCol &&d + dir[iRow][iCol] != 8&& d != 4 )
								&& !fabs(orderout[tRow][tCol]-_noData)<Eps){

								updateAtJunction(oOut,iRow,tRow,iCol,tCol,dir,orderout
											,elevout,dem,newstream,s1,s1sq,s2,s2sq,n1,n2);
							}
						}
					}
				}
				if(newstream)  // Here all paths terminated at the junction so this is a new stream
				{
					elevout[iRow][iCol]=dem[iRow][iCol];
				}
			}
			
			contribs[iRow][iCol]=-1;
	
		}
		
		if (iRow == _maxrow && iCol == _maxcol) {

			MPI_Barrier(MPI_COMM_WORLD);
			for (int i = minRow; i <= _maxrow; ++i) {
				for (int j = minCol; j <= _maxcol; ++j) {

					if (contribs[i][j] >=0) {
						bool flag=true;
						for(tRow = i - iNeighborCells; tRow <= i + iNeighborCells; tRow++) {
							for (tCol = j - iNeighborCells; tCol <= j + iNeighborCells; tCol++){				
								if(dir[tRow][tCol]>=0){
									d=dir[tRow][tCol];
									if(d1[d]+tRow ==i&& d2[d]+tCol==j &&d + dir[i][j] != 8&& d != 4 &&contribs[tRow][tCol]>-1){
										flag=false;
									}
								}
							}
						}
						if(flag){
							contribs[i][j]=0;
						}
						Termination = 0;
						num++;
					}		
				} 
            }
			//calculate t-statistic
			MPI_Barrier(MPI_COMM_WORLD);
			double gs1,gs2,gs1sq,gs2sq;
			double glen;
			long gn1,gn2;

			MPI_Allreduce(&s1,&gs1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(&s2,&gs2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(&s1sq,&gs1sq,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(&s2sq,&gs2sq,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(&n1,&gn1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(&n2,&gn2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
			MPI_Allreduce(&length,&glen,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

			double drainden=glen/(totalarea*_cellSize);
			double md1 = gs1/gn1; 
			double mdh = gs2/gn2; 
			double sd1 = sqrt((gs1sq-gn1*md1*md1)/(gn1-1)); 
			double sdh = sqrt((gs2sq-gn2*mdh*mdh)/(gn2-1)); 
			double t = (md1-mdh)/(sqrt(((gn1-1)*sd1*sd1+(gn2-1)*sdh*sdh) / (gn1+gn2-2))*sqrt(1./gn1+1./gn2));		
			opt=t;		
		}       
	}
	return true;
}




			
		
		

