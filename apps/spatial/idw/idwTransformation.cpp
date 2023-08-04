#include "idwTransformation.h"

bool IdwTransformation::Operator(const CellCoord &coord) {
    int cRow=coord.iRow();
    int cCol=coord.iCol();
    
    if( cRow == _pWorkBR->minIRow() && cCol == _pWorkBR->minICol() ){
		if( !paramInit() ){
			return false;
		}
	}
    
    const vector<SampleBlock> *sampleBlocks = _idwOperator->getSampleBlocks();
    
    CellSpace<double>& comptL = *_pComptLayer->cellSpace();

    if(_idwOperator->getMaskLayer()) {
        int maskRow= _pComptLayer->rowAtOtherLayer(_idwOperator->getMaskLayer(),cRow);
        int maskCol= _pComptLayer->colAtOtherLayer(_idwOperator->getMaskLayer(),cCol);
        int mask = (*_idwOperator->getMaskLayer()->cellSpace())[maskRow][maskCol];
        if(mask==0 | mask == _idwOperator->getMaskLayer()->metaData()->noData) {
            comptL[cRow][cCol]=_noData;
            return true;
        }
    }

    int blockRow = _idwOperator->getBlockRowIndexByCellIndex(cRow,_pComptLayer->_pMetaData->cellSize);
    int blockCol = _idwOperator->getBlockColIndexByCellIndex(cCol,_pComptLayer->_pMetaData->cellSize);
    int blockCols = _idwOperator->getBlockCols();
    int blockRows = _idwOperator->getBlockRows();

    bool isPointEnough = false;
    int searchRadius = 0; //e.g. searchRadius=1 means search the 3*3 neighborhood but without central cell
    int searchBlockNums = 0; //count of blocks need to search for this cell in comptLayer
    int nbrNums = (*sampleBlocks)[blockRow * blockCols + blockCol].samplePoints.size();
    
    while (nbrNums < _idwOperator->getNbrPoints()) //nbrPoints is the K in Eq.(7)
    {
        //search the blocks with gap cells equal to searchRadius
        for (int tRow = blockRow - searchRadius; tRow <= blockRow + searchRadius; ++tRow) {
            if (isPointEnough) break;
            if (tRow < 0 || tRow >= blockRows) continue;
            if (tRow == blockRow - searchRadius || tRow == blockRow + searchRadius) {
                for (int tCol = blockCol - searchRadius; tCol <= blockCol + searchRadius; ++tCol) {
                    if (isPointEnough) break;
                    if (tCol < 0 || tCol >= blockCols) continue;
                    // update variable nbrNums and searchBlockNums
                    nbrNums += (*sampleBlocks)[tRow * blockCols + tCol].samplePoints.size();
                    searchBlockNums++;
                    if (nbrNums >= _idwOperator->getNbrPoints())
                        isPointEnough = true;
                }
            }
        }
        if (isPointEnough) break;
        ++searchRadius;
    }
    double load_1 = searchBlockNums; //algorithms depend on specific application
    double load_2 = nbrNums * 2;
    comptL[cRow][cCol] = load_1 + load_2;
    //comptL[cRow][cCol] = 1;
    return true;
}
