#include <ogrsf_frmts.h>
#include "utility.h"
#include "reclassifyOperator.h"

void ReclassifyOperator::setInputLayer(RasterLayer<double>& inputLayer) {
    _inputLayer=&inputLayer;
    _inputLayer->newLocalNbrhood();
    Configure(_inputLayer, false);
}
bool ReclassifyOperator::Operator(const CellCoord& coord, bool operFlag) {
    int iRow = coord.iRow();
    int iCol = coord.iCol();

    CellSpace<double>& in = *_inputLayer->cellSpace();
    CellSpace<double>& out = *_outputLayer->cellSpace();
    double v = in[iRow][iCol];
    double result=0;
    if(v-_inputLayer->metaData()->noData <0.00001) {
        return true;
    }

    const double vLevelNum=10;
    double vLevels[]={0,5,30,100,200,500,1000,3500,5000,9999};

    for (int i=0; i<vLevelNum-1; i++) {
        if(v>vLevels[i]&&v<vLevels[i+1]) {
            result=i;
            break;
        }
    }

    out[iRow][iCol]=result;

    return true;
}
