#include"slopeOperator.h"


void SlopeOperator::
demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
    _pDEMNbrhood = layerD.nbrhood();
    cellSize = _pDEMLayer->_pMetaData->cellSize;
    noData = _pDEMLayer->_pMetaData->noData;
    Configure(_pDEMLayer, false);
}

void SlopeOperator::slopeLayer(RasterLayer<double>& layerD) {
    _pSlopeLayer = &layerD;
    Configure(_pSlopeLayer, false);
}

bool SlopeOperator::isTermination() {
    num--;
    if (num > 0) {
        return true;
    } else {
        return false;
    }
}

bool SlopeOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& slope = *(_pSlopeLayer->cellSpace());
    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
    //int dCellSize = _pDEMLayer->_pMetaData->cellSize;
    //int nodata = _pDEMLayer->_pMetaData->noData;

    int iRow = coord.iRow();
    int iCol = coord.iCol();

    double d[9];
    int k = 0;
    int tag = 0;

    int iRow1, iCol1;
    for (iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++) {
        for (iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++) {
            d[k] = dem[iRow1][iCol1];
            if (fabs(d[k] - noData) < Eps) {
                tag = 1;
            }
            k++;
        }
    }

    if (tag == 1) {
        slope[iRow][iCol] = noData;
        return true;
    }
    // maximum downslope
    double dm0 = (d[4] - d[0]) / sqrt(2.0);
    double dm1 = (d[4] - d[1]);
    double dm2 = (d[4] - d[2]) / sqrt(2.0);
    double dm3 = (d[4] - d[3]);
    double dm5 = (d[4] - d[5]);
    double dm6 = (d[4] - d[6]) / sqrt(2.0);
    double dm7 = (d[4] - d[7]);
    double dm8 = (d[4] - d[8]) / sqrt(2.0);
    double dmax01 = MAX(dm0, dm1);
    double dmax23 = MAX(dm2, dm3);
    double dmax56 = MAX(dm5, dm6);
    double dmax78 = MAX(dm7, dm8);
    double dmax03 = MAX(dmax01, dmax23);
    double dmax58 = MAX(dmax56, dmax78);
    double dmax = MAX(dmax03, dmax58);
    slope[iRow][iCol] = dmax / cellSize;

    return true;
}
