#include <ogrsf_frmts.h>
#include "utility.h"
#include "resampleOperator.h"

void ResampleOperator::setInputLayer(RasterLayer<double>& inputLayer) {
    _inputLayer=&inputLayer;
}
void ResampleOperator::setOutputLayer(RasterLayer<double>& outputLayer) {
    _outputLayer=&outputLayer;
    Configure(_outputLayer, false);
}
//loop output layer
bool ResampleOperator::Operator(const CellCoord& coord, bool operFlag) {
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    int myRank=GetRank();
    CellSpace<double>& in = *_inputLayer->cellSpace();
    CellSpace<double>& out = *_outputLayer->cellSpace();
    double result = 0;
    int globalRowInInputLayer = (iRow+_outputLayer->metaData()->_MBR.minIRow())*_g;
    int globalColInInputLayer = (iCol+_outputLayer->metaData()->_MBR.minICol())*_g;

    int startLocalRowInInputLayer = globalRowInInputLayer - _inputLayer->metaData()->_MBR.minIRow();
    int startLocalColInInputLayer = globalColInInputLayer - _inputLayer->metaData()->_MBR.minICol();
    
    

    //1. average
    // vector<double> cells;
    // for (int i = startLocalRowInInputLayer; i < startLocalRowInInputLayer+_g; ++i) {
    //     for (int j = startLocalColInInputLayer; j < startLocalColInInputLayer+_g; ++j) {
    //         double v = in[i][j];
    //         if(j>=_inputLayer->metaData()->_localdims.nCols()) {
    //             continue;
    //         }
    //         if(i>=_inputLayer->metaData()->_localdims.nRows()) {
    //             break;
    //         }
    //         if(!_inputLayer->isNodata(in[i][j])) {
    //             cells.push_back(in[i][j]);
    //         }
    //     }
    // }
    // if(cells.empty()) {
    //     out[iRow][iCol]=_outputLayer->_pMetaData->noData;
    //     return true;
    // }
    // double sum=0;
    // for (int i = 0; i < cells.size(); ++i) {
    //     sum+=cells[i];
    // }
    // result=sum/cells.size();

    //2. first
    if(_inputLayer->isNodata(in[startLocalRowInInputLayer][startLocalColInInputLayer])) {
        result = _outputLayer->_pMetaData->noData;
    }else {
        result = in[startLocalRowInInputLayer][startLocalColInInputLayer];
    }

    out[iRow][iCol]=result;

    return true;
}
