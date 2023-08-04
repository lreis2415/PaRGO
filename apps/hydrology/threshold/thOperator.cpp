#include"thOperator.h"

void thOperator::
scaLayer(RasterLayer<double>& layerD) {
    _pSCALayer = &layerD;
    //_pSCANbrhood = layerD.nbrhood();
    cellSize = _pSCALayer->_pMetaData->cellSize;
    noData = _pSCALayer->_pMetaData->noData;
    _xSize = _pSCALayer->_pMetaData->_localdims.nRows();
    _ySize = _pSCALayer->_pMetaData->_localdims.nCols();
    //cout<<noData<<endl;

    Configure(_pSCALayer, false);
}

void thOperator::netLayer(RasterLayer<double>& layerD) {
    _pnetLayer = &layerD;

    Configure(_pnetLayer, true);
}

bool thOperator::isTermination() {
    return flag;
}

bool thOperator::Operator(const CellCoord& coord, bool operFlag) {

    CellSpace<double>& sca = *(_pSCALayer->cellSpace());
    CellSpace<double>& net = *(_pnetLayer->cellSpace());
    int iRow = coord.iRow();
    int iCol = coord.iCol();

    //double gap = 0.0005;
	if(!fabs(sca[iRow][iCol]-noData)<Eps)
	{
		if(sca[iRow][iCol]>=th)
			net[iRow][iCol]=1;
		else
			net[iRow][iCol]=0;
	}else
		net[iRow][iCol]=noData;

    return true;
}
