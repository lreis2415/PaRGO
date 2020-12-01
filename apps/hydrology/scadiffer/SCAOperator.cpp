#include"SCAOperator.h"

#include "stdio.h"

#include<cmath>

#include<iomanip>


void SCAOperator::demLayer(RasterLayer<double>& layerD, int kc_meth2, float StepRatio2, int threadNUM2) {

    _pDEMLayer = &layerD;

    _pDEMNbrhood = layerD.nbrhood();

    kc_Meth = kc_meth2;

    StepR = StepRatio2;

    threadNum = threadNUM2;


    int xSize = _pDEMLayer->_pMetaData->_localdims.nRows(); //number of rows in this process

    int ySize = _pDEMLayer->_pMetaData->_localdims.nCols();

    zq = new double*[xSize];

    qb = new double*[xSize];

    diff = new double*[xSize];

    for (int i = 0; i < xSize; i++) {

        zq[i] = new double[ySize];

        qb[i] = new double[ySize];

        diff[i] = new double[ySize];

    }

    relDiffMin = new double[threadNum];

    relDiffMax = new double[threadNum];

    cout << "initialization is ok" << endl;

    Configure(_pDEMLayer, false);

}


void SCAOperator::SCALayer(RasterLayer<double>& layerD) {

    _pSCALayer = &layerD;

    Configure(_pSCALayer, false);

}

//void SCAOperator::tmpLayer(RasterLayer<double> &layerD) 
//
//{
//
//	_pTMPLayer = &layerD;
//
//	Configure(_pTMPLayer,false);
//
//}


bool SCAOperator::isTermination() {

    return flag;

}


bool SCAOperator::Operator(const CellCoord& coord) {

    CellSpace<double>& demL = *(_pDEMLayer->cellSpace());

    CellSpace<double>& scaL = *(_pSCALayer->cellSpace());

    //CellSpace<double> &tmpL = *(_pTMPLayer->cellSpace());

    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);

    int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;


    int xSize = _pDEMLayer->_pMetaData->_localdims.nRows(); //number of rows in this process

    int ySize = _pDEMLayer->_pMetaData->_localdims.nCols();

    int iRow = coord.iRow(); //current cell

    int iCol = coord.iCol();

    double dCellSize = _pDEMLayer->_pMetaData->cellSize; //resolution

    noData = _pDEMLayer->_pMetaData->noData;

    int mythread = omp_get_thread_num();


    if (num == 1) {

        zq[iRow][iCol] = qb[iRow][iCol] + 1.0 / 6.0 * (qb[iRow - 1][iCol] + qb[iRow + 1][iCol] + qb[iRow][iCol + 1] + qb[iRow][iCol - 1]) +

            1.0 / 36.0 * (qb[iRow - 1][iCol - 1] + qb[iRow - 1][iCol + 1] + qb[iRow + 1][iCol - 1] + qb[iRow + 1][iCol + 1]);

        diff[iRow][iCol] = zq[iRow][iCol] - demL[iRow][iCol];

        if ((iCol == ySize - 2) && ((((xSize - 2) % threadNum == 0) && (iRow == ((xSize - 2) / threadNum) * (mythread + 1))) || (((xSize - 2) % threadNum != 0) && (((mythread < (threadNum - 1)) && (iRow == ((int)(xSize - 2) / threadNum + 1) * (mythread + 1))) || ((mythread == (threadNum - 1)) && (iRow == xSize - 2)))))) {

#pragma omp barrier

            num = 2;

            Termination = 0;

#pragma omp barrier

        }

    }
    else {

        if (num == 2) {

            qb[iRow][iCol] -= 1.0 * diff[iRow][iCol] * 9.0 / 16.0;

            double tmp = fabs(demL[iRow][iCol]);

            double relDiff = tmp < 1.0 ? diff[iRow][iCol] : (diff[iRow][iCol] / tmp);

            if ((iCol == 1) && ((((xSize - 2) % threadNum == 0) && (iRow == ((xSize - 2) / threadNum) * (mythread) + 1)) || (((xSize - 2) % threadNum != 0) && (((iRow == ((int)(xSize - 2) / threadNum + 1) * mythread + 1)))))) {

                relDiffMax[mythread] = relDiff;

                relDiffMin[mythread] = relDiff;

            }

            if (relDiff < relDiffMin[mythread]) {

                relDiffMin[mythread] = relDiff;

            }
            else {

                if (relDiff > relDiffMax[mythread]) {

                    relDiffMax[mythread] = relDiff;

                }

            }

            if ((iCol == ySize - 2) && ((((xSize - 2) % threadNum == 0) && (iRow == ((xSize - 2) / threadNum) * (mythread + 1))) || (((xSize - 2) % threadNum != 0) && (((mythread < (threadNum - 1)) && (iRow == ((int)(xSize - 2) / threadNum + 1) * (mythread + 1))) || ((mythread == (threadNum - 1)) && (iRow == xSize - 2)))))) {

#pragma omp barrier

                //recalculating border

                if (mythread == 0) {

                    for (int i = 1; i < xSize - 1; i++) {

                        qb[i][0] = 2 * qb[i][1] - qb[i][2]; //left border

                        qb[i][ySize - 1] = 2 * qb[i][ySize - 2] - qb[i][ySize - 3]; //right border

                    }

                    for (int j = 1; j < ySize - 1; j++) {

                        qb[0][j] = 2 * qb[1][j] - qb[2][j]; //top border

                        qb[xSize - 1][j] = 2 * qb[xSize - 2][j] - qb[xSize - 3][j]; //bottom border

                    }

                    qb[0][0] = 4 * qb[1][1] - 2 * qb[1][2] - 2 * qb[2][1] + qb[2][2]; //top-left corner

                    qb[0][ySize - 1] = 4 * qb[1][ySize - 2] - 2 * qb[1][ySize - 3] - 2 * qb[2][ySize - 2] + qb[2][ySize - 3]; //top-right corner

                    qb[xSize - 1][0] = 4 * qb[xSize - 2][1] - 2 * qb[xSize - 2][2] - 2 * qb[xSize - 3][1] + qb[xSize - 3][2]; //bottom-left corner

                    qb[xSize - 1][ySize - 1] = 4 * qb[xSize - 2][ySize - 2] - 2 * qb[xSize - 2][ySize - 3] - 2 * qb[xSize - 3][ySize - 2] + qb[xSize - 3][ySize - 3]; //bottom-right corner

                    //whether continuing iterating approximation

                    for (int i = 1; i < threadNum; i++) {

                        if (relDiffMin[i] < relDiffMin[0]) {

                            relDiffMin[0] = relDiffMin[i];

                        }

                        if (relDiffMax[i] > relDiffMax[0]) {

                            relDiffMax[0] = relDiffMax[i];

                        }

                    }

                    if ((fabs(relDiffMax[0]) + fabs(relDiffMin[0])) >= 1e-9) {

                        num = 1;

                        Termination = 0;

                        cout << fabs(relDiffMax[0]) << " continuing iterating approximation " << fabs(relDiffMin[0]) << endl;

                        relDiffMin[0] = 0.0;

                        relDiffMax[0] = 0.0;

                    }
                    else {

                        num = 3;

                        Termination = 0;

                    }

                }

#pragma omp barrier

            }

        }
        else {

            if (num == 0) {

                qb[iRow][iCol] = 9.0 / 16.0 * demL[iRow][iCol];

                if ((iCol == ySize - 2) && ((((xSize - 2) % threadNum == 0) && (iRow == ((xSize - 2) / threadNum) * (mythread + 1))) || (((xSize - 2) % threadNum != 0) && (((mythread < (threadNum - 1)) && (iRow == ((int)(xSize - 2) / threadNum + 1) * (mythread + 1))) || ((mythread == (threadNum - 1)) && (iRow == xSize - 2)))))) {

#pragma omp barrier

                    if (mythread == 0) {

                        for (int i = 1; i < xSize - 1; i++) {

                            qb[i][0] = 2 * qb[i][1] - qb[i][2]; //left border

                            qb[i][ySize - 1] = 2 * qb[i][ySize - 2] - qb[i][ySize - 3]; //right border

                        }

                        for (int j = 1; j < ySize - 1; j++) {

                            qb[0][j] = 2 * qb[1][j] - qb[2][j]; //top border

                            qb[xSize - 1][j] = 2 * qb[xSize - 2][j] - qb[xSize - 3][j]; //bottom border

                        }

                        qb[0][0] = 4 * qb[1][1] - 2 * qb[1][2] - 2 * qb[2][1] + qb[2][2]; //top-left corner

                        qb[0][ySize - 1] = 4 * qb[1][ySize - 2] - 2 * qb[1][ySize - 3] - 2 * qb[2][ySize - 2] + qb[2][ySize - 3]; //top-right corner

                        qb[xSize - 1][0] = 4 * qb[xSize - 2][1] - 2 * qb[xSize - 2][2] - 2 * qb[xSize - 3][1] + qb[xSize - 3][2]; //bottom-left corner

                        qb[xSize - 1][ySize - 1] = 4 * qb[xSize - 2][ySize - 2] - 2 * qb[xSize - 2][ySize - 3] - 2 * qb[xSize - 3][ySize - 2] + qb[xSize - 3][ySize - 3]; //bottom-right corner

                        num = 1;

                        Termination = 0;

                    }

#pragma omp barrier

                }

            }
            else {

                if (num == 3) {

                    //check the biquadratic interpolation by comparing the value in cell center

                    //double localX = 0.5;

                    //double localY = 0.5;

                    //zq[iRow][iCol] = (4.0*localX*(1-localX)+2.0)/3.0*2.0/3.0*((qb[iRow-1][iCol]*localY*localY) + (qb[iRow][iCol]*(2.0*localY*(1-localY)+1.0)) + (qb[iRow+1][iCol]*(1-localY)*(1-localY))) +

                    //	2.0*localX*localX/3.0*2.0/3.0*((qb[iRow-1][iCol+1]*localY*localY) + (qb[iRow][iCol+1]*(2.0*localY*(1-localY)+1.0)) + (qb[iRow+1][iCol+1]*(1-localY)*(1-localY))) +

                    //	2.0*(1.0-localX)*(1.0-localX)/3.0*2.0/3.0*((qb[iRow-1][iCol-1]*localY*localY) + (qb[iRow][iCol-1]*(2.0*localY*(1-localY)+1.0)) + (qb[iRow+1][iCol-1]*(1-localY)*(1-localY)));

                    //zq[iRow][iCol] = q2(localX)*(qb[iRow-1][iCol]*q1(localY)+qb[iRow][iCol]*q2(localY)+qb[iRow+1][iCol]*q3(localY)) +

                    //				q1(localX)*(qb[iRow-1][iCol+1]*q1(localY)+qb[iRow][iCol+1]*q2(localY)+qb[iRow+1][iCol+1]*q3(localY)) +

                    //				q3(localX)*(qb[iRow-1][iCol-1]*q1(localY)+qb[iRow][iCol-1]*q2(localY)+qb[iRow+1][iCol-1]*q3(localY));

                    vector<double> planCur;

                    vector<double> stepL;

                    int curRow = iRow;

                    int curCol = iCol; //which cell the current point is in

                    double curX = 0.5 * dCellSize;

                    double curY = 0.5 * dCellSize; //the local coordinate of current point

                    double prevAsp = 0.0;

                    double nextStepL = StepR * dCellSize;

                    double curStepL = StepR * dCellSize;

                    bool upFlag = true;

                    int expNum = 0;

                    while (upFlag) {

                        //double globalx = -400+dCellSize*curCol+curX;

                        //double globaly = -400+dCellSize*(xSize-curRow-1)+curY;

                        //double globalx = 437168.40625+(double)(dCellSize*curCol)+curX;

                        //double globaly = 5417476.5+(double)dCellSize*(xSize-curRow-1)+curY;

                        //if( (iRow==450)&&(iCol==500) )

                        //if ((iRow==690)&&(iCol==460))

                        //cout<<"globalx "<<globalx<<" globaly "<<globaly<<" tmpL[curRow][curCol] "<<tmpL[curRow][curCol]<<endl;

                        //{printf("%3f,%3f",globalx,globaly);

                        //cout<<endl;}

                        //cout<<"tmpL[curRow][curCol] "<<tmpL[curRow][curCol]<<endl;

                        if ((curRow < 1) || (curRow > xSize - 2) || (curCol < 1) || (curCol > ySize - 2)) {

                            upFlag = false;

                        }
                        else {

                            //whether noData cells exists in its neighborhood

                            for (int iRow1 = curRow - iNeighborCells; iRow1 <= curRow + iNeighborCells; iRow1++) {

                                if (!upFlag) {

                                    break;

                                }

                                for (int iCol1 = curCol - iNeighborCells; iCol1 <= curCol + iNeighborCells; iCol1++) {

                                    if (fabs(demL[iRow1][iCol1] - noData) < Eps) {

                                        upFlag = false;

                                        break;

                                    }

                                }

                            }

                            //hilltop

                            if (upFlag) {

                                upFlag = false;

                                for (int iRow1 = curRow - iNeighborCells; iRow1 <= curRow + iNeighborCells; iRow1++) {

                                    if (upFlag) {

                                        break;

                                    }

                                    for (int iCol1 = curCol - iNeighborCells; iCol1 <= curCol + iNeighborCells; iCol1++) {

                                        if (demL[iRow1][iCol1] - demL[curRow][curCol] > Eps) {

                                            upFlag = true;

                                            break;

                                        }

                                    }

                                }

                                if (!upFlag) {

                                    scaL[curRow][curCol] = 0.0;

                                    //calculating when the endpoint is hilltop

                                    if (stepL.size() > 0) {

                                        double z[9];

                                        int dir = 0;

                                        for (int iRow1 = curRow - iNeighborCells; iRow1 <= curRow + iNeighborCells; iRow1++) {

                                            for (int iCol1 = curCol - iNeighborCells; iCol1 <= curCol + iNeighborCells; iCol1++) {

                                                z[dir] = qb[iRow1][iCol1];

                                                dir++;

                                            }

                                        }

                                        double p, q, r, s, t;

                                        double percurX = curX / dCellSize;

                                        double percurY = curY / dCellSize;

                                        p = 4.0 * (1 - 2.0 * percurX) / 3.0 * (z[1] * q1(percurY) + z[4] * q2(percurY) + z[7] * q3(percurY)) +

                                            4.0 * percurX / 3.0 * (z[2] * q1(percurY) + z[5] * q2(percurY) + z[8] * q3(percurY)) +

                                            4.0 * (percurX - 1.0) / 3.0 * (z[0] * q1(percurY) + z[3] * q2(percurY) + z[6] * q3(percurY));

                                        r = -8.0 / 3.0 * (z[1] * q1(percurY) + z[4] * q2(percurY) + z[7] * q3(percurY)) +

                                            4.0 / 3.0 * (z[2] * q1(percurY) + z[5] * q2(percurY) + z[8] * q3(percurY)) +

                                            4.0 / 3.0 * (z[0] * q1(percurY) + z[3] * q2(percurY) + z[6] * q3(percurY));

                                        q = 4.0 * (1 - 2.0 * percurY) / 3.0 * (z[5] * q1(percurX) + z[4] * q2(percurX) + z[3] * q3(percurX)) +

                                            4.0 * percurY / 3.0 * (z[2] * q1(percurX) + z[1] * q2(percurX) + z[0] * q3(percurX)) +

                                            4.0 * (percurY - 1.0) / 3.0 * (z[8] * q1(percurX) + z[7] * q2(percurX) + z[6] * q3(percurX));

                                        t = -8.0 / 3.0 * (z[5] * q1(percurX) + z[4] * q2(percurX) + z[3] * q3(percurX)) +

                                            4.0 / 3.0 * (z[2] * q1(percurX) + z[1] * q2(percurX) + z[0] * q3(percurX)) +

                                            4.0 / 3.0 * (z[8] * q1(percurX) + z[7] * q2(percurX) + z[6] * q3(percurX));

                                        s = 4.0 * (1 - 2.0 * percurY) / 3.0 * (4.0 * (1 - 2.0 * percurX) / 3.0 * z[4] + 4.0 * percurX / 3.0 * z[5] + 4.0 * (percurX - 1.0) / 3.0 * z[3]) +

                                            4.0 * percurY / 3.0 * (4.0 * (1 - 2.0 * percurX) / 3.0 * z[1] + 4.0 * percurX / 3.0 * z[2] + 4.0 * (percurX - 1.0) / 3.0 * z[0]) +

                                            4.0 * (percurY - 1.0) / 3.0 * (4.0 * (1 - 2.0 * percurX) / 3.0 * z[7] + 4.0 * percurX / 3.0 * z[8] + 4.0 * (percurX - 1.0) / 3.0 * z[6]);

                                        double tmpCur = -((q * q * r - 2 * p * q * s + p * p * t) / (sqrt((p * p + q * q) * (p * p + q * q) * (p * p + q * q)))) / dCellSize;

                                        planCur.push_back(tmpCur);

                                        stepL.push_back(curStepL);

                                    }

                                }

                            }

                        }

                        if (upFlag) {

                            double z[9];

                            int dir = 0;

                            for (int iRow1 = curRow - iNeighborCells; iRow1 <= curRow + iNeighborCells; iRow1++) {

                                for (int iCol1 = curCol - iNeighborCells; iCol1 <= curCol + iNeighborCells; iCol1++) {

                                    z[dir] = qb[iRow1][iCol1];

                                    dir++;

                                }

                            }

                            double p, q, r, s, t;

                            double percurX = curX / dCellSize;

                            double percurY = curY / dCellSize;

                            p = 4.0 * (1 - 2.0 * percurX) / 3.0 * (z[1] * q1(percurY) + z[4] * q2(percurY) + z[7] * q3(percurY)) +

                                4.0 * percurX / 3.0 * (z[2] * q1(percurY) + z[5] * q2(percurY) + z[8] * q3(percurY)) +

                                4.0 * (percurX - 1.0) / 3.0 * (z[0] * q1(percurY) + z[3] * q2(percurY) + z[6] * q3(percurY));

                            r = -8.0 / 3.0 * (z[1] * q1(percurY) + z[4] * q2(percurY) + z[7] * q3(percurY)) +

                                4.0 / 3.0 * (z[2] * q1(percurY) + z[5] * q2(percurY) + z[8] * q3(percurY)) +

                                4.0 / 3.0 * (z[0] * q1(percurY) + z[3] * q2(percurY) + z[6] * q3(percurY));

                            q = 4.0 * (1 - 2.0 * percurY) / 3.0 * (z[5] * q1(percurX) + z[4] * q2(percurX) + z[3] * q3(percurX)) +

                                4.0 * percurY / 3.0 * (z[2] * q1(percurX) + z[1] * q2(percurX) + z[0] * q3(percurX)) +

                                4.0 * (percurY - 1.0) / 3.0 * (z[8] * q1(percurX) + z[7] * q2(percurX) + z[6] * q3(percurX));

                            t = -8.0 / 3.0 * (z[5] * q1(percurX) + z[4] * q2(percurX) + z[3] * q3(percurX)) +

                                4.0 / 3.0 * (z[2] * q1(percurX) + z[1] * q2(percurX) + z[0] * q3(percurX)) +

                                4.0 / 3.0 * (z[8] * q1(percurX) + z[7] * q2(percurX) + z[6] * q3(percurX));

                            s = 4.0 * (1 - 2.0 * percurY) / 3.0 * (4.0 * (1 - 2.0 * percurX) / 3.0 * z[4] + 4.0 * percurX / 3.0 * z[5] + 4.0 * (percurX - 1.0) / 3.0 * z[3]) +

                                4.0 * percurY / 3.0 * (4.0 * (1 - 2.0 * percurX) / 3.0 * z[1] + 4.0 * percurX / 3.0 * z[2] + 4.0 * (percurX - 1.0) / 3.0 * z[0]) +

                                4.0 * (percurY - 1.0) / 3.0 * (4.0 * (1 - 2.0 * percurX) / 3.0 * z[7] + 4.0 * percurX / 3.0 * z[8] + 4.0 * (percurX - 1.0) / 3.0 * z[6]);

                            double tmpCur = 1.0;

                            tmpCur = -((q * q * r - 2 * p * q * s + p * p * t) / (sqrt((p * p + q * q) * (p * p + q * q) * (p * p + q * q)))) / dCellSize;

                            planCur.push_back(tmpCur);

                            if (!(stepL.empty())) {

                                //calculating aspect,north is zero degree in a clockwise direction

                                double tmpAsp = 0.0;

                                if (fabs(p - 0) < Eps) {

                                    if (fabs(q - 0) < Eps) {

                                        //Cells in the input grid of zero slope (flat) are calculated by D8

                                        double tmpDiffer = 0.0;

                                        int tmpLowNum = 4;

                                        for (int i = 0; i < 9; i++) {

                                            if ((i % 2 == 0) && (i != 4) && ((z[4] - z[i]) / sqrt(2.0) > tmpDiffer)) {

                                                tmpDiffer = (z[4] - z[i]) / sqrt(2.0);

                                                tmpLowNum = i;

                                            }
                                            else {

                                                if ((i % 2 != 0) && ((z[4] - z[i]) > tmpDiffer)) {

                                                    tmpDiffer = z[4] - z[i];

                                                    tmpLowNum = i;

                                                }

                                            }

                                        }

                                        switch (tmpLowNum) {

                                        case 0: tmpAsp = 315;
                                            break;

                                        case 1: tmpAsp = 0;
                                            break;

                                        case 2: tmpAsp = 45;
                                            break;

                                        case 3: tmpAsp = 270;
                                            break;

                                        case 4: tmpAsp = false;
                                            break;

                                            //case 4:cout<<"start mid-point central aspect is wrong for it's a pit for iRow="<<iRow<<" iCol="<<iCol<<endl;break;

                                        case 5: tmpAsp = 90;
                                            break;

                                        case 6: tmpAsp = 225;
                                            break;

                                        case 7: tmpAsp = 180;
                                            break;

                                        case 8: tmpAsp = 135;
                                            break;

                                        default: cout << "start mid-point central aspect is wrong(impossible)" << endl;
                                            break;

                                        }

                                        //cout<<"p && q is nearly zero"<<endl;

                                    }
                                    else {

                                        tmpAsp = 90.0 + 90.0 * q / fabs(q);

                                        //cout<<"p is nearly zero"<<endl;

                                    }

                                }
                                else {

                                    tmpAsp = 180.0 - 57.29578 * atan(q / p) + 90.0 * p / fabs(p);

                                }

                                //compare the aspect to find next integration step length

                                //if(fabs(tmpAsp-prevAsp) <= 0.8){

                                //	nextStepL = 1.5*curStepL;

                                //	if(nextStepL/dCellSize > 10){

                                //		nextStepL = dCellSize * 10;	//防止步长过大

                                //	}

                                //}else{

                                //	nextStepL = 0.5*curStepL;

                                //	if(nextStepL/dCellSize < 0.1){

                                //		nextStepL = dCellSize * 0.1;	//防止步长过小

                                //		cout<<"too short step length for "<<iRow<<" "<<iCol<<" with "<<stepL.size()<<endl;

                                //	}

                                //}

                                stepL.push_back(nextStepL);

                                //specially processing for infinite loop

                                if ((fabs(tmpAsp - prevAsp) > 179) && (fabs(tmpAsp - prevAsp) < 181)) {

                                    expNum++;

                                }
                                else {

                                    expNum = 0;

                                }

                                if (expNum > 10) {

                                    upFlag = false;

                                }

                                prevAsp = tmpAsp;

                            }
                            else {

                                stepL.push_back(curStepL);

                                //calculating aspect for start point,north is zero degree in a clockwise direction

                                if (fabs(p - 0.0) < Eps) {

                                    if (fabs(q - 0.0) < Eps) {

                                        //Cells in the input grid of zero slope (flat) are approximately calculated by D8

                                        double tmpDiffer = 0.0;

                                        int tmpLowNum = 4;

                                        for (int i = 0; i < 9; i++) {

                                            if ((i % 2 == 0) && (i != 4) && ((z[4] - z[i]) / sqrt(2.0) > tmpDiffer)) {

                                                tmpDiffer = (z[4] - z[i]) / sqrt(2.0);

                                                tmpLowNum = i;

                                            }
                                            else {

                                                if ((i % 2 != 0) && ((z[4] - z[i]) > tmpDiffer)) {

                                                    tmpDiffer = z[4] - z[i];

                                                    tmpLowNum = i;

                                                }

                                            }

                                        }

                                        switch (tmpLowNum) {

                                        case 0: prevAsp = 315;
                                            break;

                                        case 1: prevAsp = 0;
                                            break;

                                        case 2: prevAsp = 45;
                                            break;

                                        case 3: prevAsp = 270;
                                            break;

                                        case 4: upFlag = false;
                                            break;

                                            //case 4:cout<<"start mid-point central aspect is wrong for it's a pit for iRow="<<iRow<<" iCol="<<iCol<<endl;break;

                                        case 5: prevAsp = 90;
                                            break;

                                        case 6: prevAsp = 225;
                                            break;

                                        case 7: prevAsp = 180;
                                            break;

                                        case 8: prevAsp = 135;
                                            break;

                                        default: cout << "start mid-point central aspect is wrong(impossible)" << endl;
                                            break;

                                        }

                                        //cout<<"p && q is nearly zero"<<endl;

                                    }
                                    else {

                                        prevAsp = 90.0 + 90.0 * q / fabs(q);

                                        //cout<<"p is nearly zero"<<endl;

                                    }

                                }
                                else {

                                    prevAsp = 180.0 - 57.29578 * atan(q / p) + 90.0 * p / fabs(p);

                                }

                            }

                            //求当前点的上游点所在栅格及其局部坐标，步长curStepL

                            int upRow, upCol;

                            double upx, upy;

                            upx = curX - curStepL * sin(prevAsp / 180.0 * 3.1415926);

                            upy = curY - curStepL * cos(prevAsp / 180.0 * 3.1415926);

                            if (upx < 0) {

                                upCol = curCol - (int)(fabs(upx) / dCellSize) - 1;

                                if (upy < 0) {

                                    upRow = curRow + (int)(fabs(upy) / dCellSize) + 1;

                                }
                                else {

                                    upRow = curRow - (int)(upy / dCellSize);

                                }

                            }
                            else {

                                upCol = curCol + (int)(upx / dCellSize);

                                if (upy < 0) {

                                    upRow = curRow + (int)(fabs(upy) / dCellSize) + 1;

                                }
                                else {

                                    upRow = curRow - (int)(upy / dCellSize);

                                }

                            }

                            upx = upx - (upCol - curCol) * dCellSize;

                            upy = upy - (curRow - upRow) * dCellSize;


                            curStepL = nextStepL;

                            curRow = upRow;

                            curCol = upCol;

                            curX = upx;

                            curY = upy;

                        }

                    }

                    //integral downslope

                    double a0, a, k0, k1, k, l, c0, c1;

                    a0 = 0.0;

                    int segNum = stepL.size();

                    if (segNum > 0) {

                        for (int i = segNum - 1; i > 0; i--) {

                            int tmpMeth = kc_Meth;

                            l = stepL[i - 1];

                            k0 = planCur[i];

                            k1 = planCur[i - 1];

                            if (tmpMeth == 2) {

                                if ((fabs(k0 - 0) > Eps) && (fabs(k1 - 0) > Eps)) {

                                    c0 = 1.0 / k0;

                                    c1 = (1.0 / k1 - 1.0 / k0) / l;

                                    if (fabs(1 + c1) < 1e-9 || fabs(c1) < 1e-9 || k0 * k1 < 0) {

                                        tmpMeth = 1;

                                    }
                                    else {

                                        a = (c0 + c1 * l) / (1 + c1) + (a0 - c0 / (1 + c1)) * pow((c0 + c1 * l) / c0, -1 / c1);

                                    }

                                }
                                else {

                                    tmpMeth = 1;

                                }

                            }

                            if (tmpMeth == 1) {

                                k = (k0 + k1) / 2.0;

                                if ((k * l) > 10 * dCellSize) {
                                    //special case that (k*l) is too big

                                    a = 1.0 / k;

                                    //cout<<"(k*l) is too big, so 1.0/k approximation, k*l = "<<k*l<<" a = "<<a<<endl;

                                }
                                else {

                                    if (fabs(k * l) < 0.0001) {

                                        //cout<<"fabs(k*l) is too small, so linear approximation, k*l = "<<k*l<<" a = "<<a<<endl;

                                        a = a0 + l;

                                    }
                                    else {

                                        a = 1. / k - (1. / k - a0) * exp(-k * l);

                                    }

                                }

                            }


                            if (a > 30000) {
                                //the extreme SCA value should be determined by practical application

                                a = noData;

                                i = 0;

                                //cout<<a<<"too large, no further calculation"<<iRow<<","<<iCol<<endl;

                            }

                            a0 = a;

                            //if ((iRow==690)&&(iCol==460))

                            //cout<<a0<<endl;

                        }

                        scaL[iRow][iCol] = a;

                        vector<double>().swap(planCur);

                        vector<double>().swap(stepL);


                    }

                }
                else {

                    cout << "something is wrong with iteration" << endl;

                }

            }

        }

    }


    return true;

}
