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
    return num > 0;
    // The if-else can be simplified by `return num > 0` --LJ
    //if(num > 0)
    //{
    //	return true;
    //}
    //else
    //{
    //	return false;
    //}
}


bool MFDOperator::Operator(const CellCoord &coord,bool operFlag)
{	
    CellSpace<double> &demL = *(_pDEMLayer->cellSpace());
    vector<CellSpace<double>* > weightLs(8);
    for( int i=0; i<8; ++i ){
        weightLs[i] = _weightLayers[i]->cellSpace();
    }
    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
	
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    // TODO: Considering should floor() or ceil() be necessary? I think this should be one function of `Neighborhood` class. --LJ
    //       Or `int GetRowNum() { return sqrt((double)(_pDEMNbrhood->size() - 1) / 2); }`? --line 29 in LEOperator.h
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	
    // TODO: Is it possible to add a function to checkout NODATA for `CellSpace` like demL.isNodata(iRow, iCol)? --LJ
    if (fabs(demL[iRow][iCol] - _noData) <= Eps) {
        for (int i = 0; i < 8; ++i) {
            (*weightLs[i])[iRow][iCol] = _noData;
        }
        return true;
    }
	
    //if(fabs(demL[iRow][iCol] - _noData) > Eps){
        for( int i = 0; i < 8; ++i ){
            (*weightLs[i])[iRow][iCol] = 0.0;
        }
		//double d[9];//邻域栅格高程
		int dir = 0;//邻域栅格编号
		double maxDiffVal = 0.0;//最大坡降正切值
		double pe = 0.0;//公式中的坡度指数p
        double l[9];             //weighted factor of contour line; parameter in equation
        double diffVal[9];

        // TODO: Some commonly used constants should be predefined in a high precision, e.g., SQ2 = 1.4142135623730951
        // Refer to the origin implementation in SimDTA/src/modMFD.bas, start from line 457. --LJ
        for(int tRow = iRow-iNeighborCells; tRow <= iRow+iNeighborCells; tRow++){
            for(int tCol = iCol-iNeighborCells; tCol <= iCol+iNeighborCells; tCol++){
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
                            diffVal[dir] = (demL[iRow][iCol]-demL[tRow][tCol])/(((double)sqrt((double)2)) * _cellSize);//1.41->sqrt(2) --FXC
                            l[dir] = ((double)sqrt((double)2)) / 4 * _cellSize; //? -- sqrt(2)/4  Should be precisely predefined! --LJ
                        }else{
                            diffVal[dir] = (demL[iRow][iCol]-demL[tRow][tCol])/_cellSize;
                            l[dir] = 0.5 * _cellSize;
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
		pe = (10 - slpExp) * maxDiffVal + slpExp; //1->10 --FXC    		
		double sum = 0.0;
        for(int i = 0;i < 9;i++){
            sum += pow(diffVal[i], pe) * l[i];
        }
        if(sum > 0){
            for(int i = 0;i < 4;i++){
                (*weightLs[i])[iRow][iCol] = (pow(diffVal[i], pe) * l[i]) / sum ;
            }
            for(int i = 5;i < 9;i++){
                (*weightLs[i-1])[iRow][iCol] = (pow(diffVal[i], pe) * l[i]) / sum ;
            }
        }
    //}
    //else{ // This else should be excluded at the front, which can save a indent of the above code. --LJ
    //    for( int i=0; i<8; ++i ){
    //        (*weightLs[i])[iRow][iCol] = _noData;
    //    }
    //}

    return true;
}
