#ifndef RUNOPERATOR_H
#define RUNOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include "transformation.h"
#include "utility.h"
#include <cmath>
#include <functional>

using namespace GPRO;

#define Eps 0.0000001

class FCMOperator : public RasterOperator<double> {
public:
    FCMOperator()
        : RasterOperator<double>(),
          _iterNum(0), flag(true),
          block(0), nodataNums(0), dFlag(true), subval(0.0), totval(0.0), oldtval(0.0), partitionCoef(0.0), entropy(0), totpartitionCoef(0), totentropy(0), reduceTime(0) {
    }

    ~FCMOperator();

    void initialization(int iNum, int cNum, int maxIter, double toler, double m);
    void inputLayer(vector<RasterLayer<double>*> layerD);
    void fcmLayer(RasterLayer<double>& layerD);
    void degLayer(vector<RasterLayer<double>*> layerD);

    virtual bool isTermination();

    bool isNoData(double value, int rasterLayerIndex);
    bool allNoDataAt(int row, int col);
    bool anyNoDataAt(int row, int col);
    inline bool firstNoDataAt(int row, int col);
    void createRandomIdx(int nums, int range, int* randomIdx);
    void fnDistance(int curRow, int curCol, double* pInputVal);
    void InitDegree(int curRow, int curCol);
    void initRandomClusterCenters(double* clusterCenters);
    void assignMaxMembershipDegrees();
    virtual bool Operator(const CellCoord& coord, bool operFlag);

    double reduceTime;
private:
    int _cellSize;
    int _xSize, _ySize; ///<local size
    int _nRows, _nCols; ///<global size
    double _noData;
    int _rank;
    int _iterNum; ///<iteration number
    bool flag;
    double preExpTime;
    double tmpSumTime1, tmpSumTime2;
protected:
    int clusterNum;
    int maxIteration;
    double tolerance; ///<threshold of iterations
    double weight;
    int imageNum; ///<input layer num

    int block;
    int nodataNums;
    bool dFlag;
    double* centerVal;
    int* centerIndex;
    double*** degree;
    double*** dist;
    double subval;
    double totval;
    double oldtval;
    double* sumNumerator;
    double* sumDenominator;
    double* totNumerator;
    double* totDenominator;
    double partitionCoef;
    double entropy;
    double totpartitionCoef; ///<total
    double totentropy; ///<total

    vector<RasterLayer<double>*> _vInputLayer;
    RasterLayer<double>* _pFCMLayer;
    vector<RasterLayer<double>*> _vDegLayer;

    Neighborhood<double>* _pDEMNbrhood;

};

#endif
