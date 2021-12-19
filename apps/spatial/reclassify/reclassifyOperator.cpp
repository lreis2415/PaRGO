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


    for (int i=0; i<levels->size()-1; i++) {
        if(v>=(*levels)[i]&&v<(*levels)[i+1]) {
            result=i;
            break;
        }
    }

    out[iRow][iCol]=result;

    return true;
}

