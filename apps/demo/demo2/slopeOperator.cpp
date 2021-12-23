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
    return num > 0;
}

bool SlopeOperator::Operator(const CellCoord& coord, bool operFlag) {
    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& slope = *(_pSlopeLayer->cellSpace());
    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;

    int iRow = coord.iRow();
    int iCol = coord.iCol();

    double d[9];
    int k = 0;
    int tag = 0;

    for (int iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++) {
        for (int iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++) {
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

    double dx = (d[4] - d[3]) / cellSize;
    double dy = (d[4] - d[7]) / cellSize;
    slope[iRow][iCol] = sqrt(dx * dx + dy * dy);

    return true;
}
