#include"mfdOperator.h"
//#include "math.h"
#define Eps 0.000001

void MFDOperator::
demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  _cellSize = _pDEMLayer->_pMetaData->cellSize;
  _noData = _pDEMLayer->_pMetaData->noData;
  Configure(_pDEMLayer, false);
}

void MFDOperator::weightLayers(vector<RasterLayer<double> *> &layers) 
{
  for( int i=0; i<8; ++i ){
	  _weightLayers.push_back(layers[i]);
	  Configure(_weightLayers[i],false);
  }
}

bool MFDOperator::isTermination()
{
	num--;
	if(num > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool MFDOperator::Operator(const CellCoord &coord,bool operFlag)
{
	CellSpace<double> &demL = *(_pDEMLayer->cellSpace());
	vector<CellSpace<double>* > weightLs;
	for( int i=0; i<8; ++i ){
		weightLs.push_back(_weightLayers[i]->cellSpace());
	}
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
	
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	

	if(fabs(demL[iRow][iCol] - _noData) > Eps){
		for( int i=0; i<8; ++i ){
			(*weightLs[i])[iRow][iCol] = 0.0;
		}
		//cout<<rowNum<<colNum<<endl;
		//double d[9];//邻域栅格高程
		int dir = 0;//邻域栅格编号
		double maxDiffVal = 0.0;//最大坡降正切值
		double pe = 0.0;//公式中的坡度指数p
		double l[9];//weighted factor of contour line; parameter in equation
		double diffVal[9];

		for(int tRow=iRow-iNeighborCells; tRow<=iRow+iNeighborCells; tRow++){
			for(int tCol=iCol-iNeighborCells; tCol<=iCol+iNeighborCells; tCol++){
				if(fabs(demL[tRow][tCol] - _noData) <= Eps){
					//if(dir<4)
					//	*(weightLs[dir])[iRow][iCol] = 0;
					//else
					//	*(weightLs[dir-1])[iRow][iCol] = 0;
					diffVal[dir] = 0.0;
					l[dir] = 0.0;
				}else{
					if((demL[iRow][iCol] <= demL[tRow][tCol])){
						diffVal[dir] = 0.0;
						l[dir] = 0.0;
					}else{
						if(dir==0||dir==2||dir==6||dir==8){
							diffVal[dir] = (demL[iRow][iCol]-demL[tRow][tCol])/(1.414*_cellSize);
							l[dir] = 0.354*_cellSize;//?
						}else{
							diffVal[dir] = (demL[iRow][iCol]-demL[tRow][tCol])/_cellSize;
							l[dir] = 0.5*_cellSize;
						}
						if (diffVal[dir] > maxDiffVal)	
							maxDiffVal = diffVal[dir];
					}
				}
				dir++;
			}
		}
		if(maxDiffVal > 1)
			maxDiffVal = 1;
		pe = (1-slpExp)*maxDiffVal + slpExp;
		double sum = 0.0;
		for(int i=0;i<9;i++){
			sum += pow(diffVal[i],pe)*l[i];
		}
		if(sum > 0){
			for(int i=0;i<4;i++){
				(*weightLs[i])[iRow][iCol] = (pow(diffVal[i],pe)*l[i])/sum ;
			}
			for(int i=5;i<9;i++){
				(*weightLs[i-1])[iRow][iCol] = (pow(diffVal[i],pe)*l[i])/sum ;
			}
		}
	}else{
		for( int i=0; i<8; ++i ){
			(*weightLs[i])[iRow][iCol] = _noData;
		}
	}

	return true;
}

