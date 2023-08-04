#include"pitRemoveOperator.h"

void PitRemoveOperator::
demLayer(RasterLayer<double>& layerD) {
    _pDEMLayer = &layerD;
    _pDEMNbrhood = layerD.nbrhood();
    cellSize = _pDEMLayer->_pMetaData->cellSize;
    noData = _pDEMLayer->_pMetaData->noData;
    _xSize = _pDEMLayer->_pMetaData->_localdims.nRows();
    _ySize = _pDEMLayer->_pMetaData->_localdims.nCols();
    //cout<<noData<<endl;

    Configure(_pDEMLayer, false);
}

void PitRemoveOperator::wdemLayer(RasterLayer<double>& layerD) {
    _pwDEMLayer = &layerD;

    Configure(_pwDEMLayer, true);
}

bool PitRemoveOperator::isTermination() {
    return flag;
}

bool PitRemoveOperator::Operator(const CellCoord& coord, bool operFlag) {

    CellSpace<double>& dem = *(_pDEMLayer->cellSpace());
    CellSpace<double>& wdem = *(_pwDEMLayer->cellSpace());
    int iRow = coord.iRow();
    int iCol = coord.iCol();

    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
    int iNeighborCells = (int)sqrt((double)nbrhoodD.size()) / 2;

    double gap1 = 0.0005;
	double gap2 = gap1 * ((double)sqrt((double)2));

    int i, j;
    if (num == 0) {
        if (fabs(dem[iRow][iCol] - noData) > Eps) {
            wdem[iRow][iCol] = 10000.0;
        }
        else {
            wdem[iRow][iCol] = noData;
        }
        if ((iRow == _xSize - 2) && (iCol == _ySize - 2)) {
            MPI_Barrier(MPI_COMM_WORLD);
            num = 1;
            Termination = 0;
        }
    }
    else {
        if (fabs(dem[iRow][iCol] - noData) <= Eps || wdem[iRow][iCol] <= dem[iRow][iCol]) {
            return true;
        }


		int dir = 0;
		double gap = 0.;
        for (i = iRow - iNeighborCells; i <= iRow + iNeighborCells; i++) {
            for (j = iCol - iNeighborCells; j <= iCol + iNeighborCells; j++) {
				if (dir % 2 == 0) // different gap in different directions --FXC
					gap = gap2;
				else
					gap = gap1;
                if ((dem[iRow][iCol] >= (wdem[i][j] + gap)) || fabs(wdem[i][j] - noData) < Eps) {
                    wdem[iRow][iCol] = dem[iRow][iCol];
                    Termination = 0;
                }else if (wdem[iRow][iCol] > (wdem[i][j] + gap)) {
                    wdem[iRow][iCol] = wdem[i][j] + gap;
                    Termination = 0;
                }
				dir++;
            }
        }
    }

    return true;
}
