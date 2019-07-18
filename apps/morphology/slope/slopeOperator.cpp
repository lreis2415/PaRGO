#include"slopeOperator.h"

SlopeAlgor GetSlopeAlgorithm(const string& arg) {
    if (StringMatch(arg, "FD")) {
        return FD;
    }
    if (StringMatch(arg, "FFD")) {
        return FFD;
    }
    if (StringMatch(arg, "MD")) {
        return MD;
    }
    if (StringMatch(arg, "SD")) {
        return SD;
    }
    if (StringMatch(arg, "SFD")) {
        return SFD;
    }
    if (StringMatch(arg, "TFD")) {
        return TFD;
    }
    if (StringMatch(arg, "TFDW")) {
        return TFDW;
    }
    return FD;
}

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

void SlopeOperator::calcAlgorithm(SlopeAlgor algor) {
    calcAlgor = algor;
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

    // What if iNeighborCells greater than 1? then k will greater than 8! By ljzhu
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

    double dx;
    double dy;

    if (calcAlgor == FD) { // Third-order finite difference weighted by reciprocal of squared distance
        dx = (d[8] + 2 * d[5] + d[2] - d[6] - 2 * d[3] - d[0]) / (8.0 * cellSize);
        dy = (d[2] + 2 * d[1] + d[0] - d[6] - 2 * d[7] - d[8]) / (8.0 * cellSize);
        slope[iRow][iCol] = sqrt(dx * dx + dy * dy);
    } else if (calcAlgor == FFD) { //Frame finite difference
        dx = (d[8] - d[6] + d[2] - d[0]) / (4.0 * cellSize);
        dy = (d[0] - d[6] + d[2] - d[8]) / (4.0 * cellSize);
        slope[iRow][iCol] = sqrt(dx * dx + dy * dy);
    } else if (calcAlgor == MD) { // Maximum downslope
        double max_diff = 0.;
        double hori_dist = cellSize;
        for (int tmpr = iRow - iNeighborCells; tmpr <= iRow + iNeighborCells; tmpr++) {
            for (int tmpc = iCol - iNeighborCells; tmpc <= iCol + iNeighborCells; tmpc++) {
                if (tmpr == iRow && tmpc == iCol) continue;
                if (max_diff < d[4] - dem[tmpr][tmpc]) {
                    max_diff = d[4] - dem[tmpr][tmpc];
                }
                hori_dist = cellSize * sqrt(double((tmpr - iRow) * (tmpr - iRow) + (tmpc - iCol) * (tmpc - iCol)));
            }
        }
        slope[iRow][iCol] = max_diff / hori_dist;

        /*double dm0 = (d[4] - d[0]) / sqrt(2.0);
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
        slope[iRow][iCol] = dmax / cellSize;*/
    } else if (calcAlgor == SD) { // Simple difference
        dx = (d[4] - d[3]) / cellSize;
        dy = (d[4] - d[7]) / cellSize;
        slope[iRow][iCol] = sqrt(dx * dx + dy * dy);
    } else if (calcAlgor == SFD) { // Second-order finite difference
        dx = (d[5] - d[3]) / (2.0 * cellSize);
        dy = (d[1] - d[7]) / (2.0 * cellSize);
        slope[iRow][iCol] = sqrt(dx * dx + dy * dy);
    } else if (calcAlgor == TFD) { // Third-order finite difference
        dx = (d[8] + d[5] + d[2] - d[0] - d[3] - d[6]) / (6.0 * cellSize);
        dy = (d[2] + d[1] + d[0] - d[6] - d[7] - d[6]) / (6.0 * cellSize);
        slope[iRow][iCol] = sqrt(dx * dx + dy * dy);
    } else if (calcAlgor == TFDW) { // Third-order finite difference weighted by reciprocal of distance)
        dx = (d[8] + sqrt(2.0) * d[5] + d[2] - d[6] - sqrt(2.0) * d[3] - d[0]) / ((4.0 + sqrt(8.0)) * cellSize);
        dy = (d[2] + sqrt(2.0) * d[1] + d[0] - d[6] - sqrt(2.0) * d[7] - d[8]) / ((4.0 + sqrt(8.0)) * cellSize);
        slope[iRow][iCol] = sqrt(dx * dx + dy * dy);
    }

    return true;
}
