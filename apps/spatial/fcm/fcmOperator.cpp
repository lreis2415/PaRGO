#include "fcmOperator.h"

FCMOperator::~FCMOperator() {
    delete centerVal;
    delete centerIndex;
    delete sumNumerator;
    delete sumDenominator;
    delete totNumerator;
    delete totDenominator;
    for (int i = 0; i < clusterNum; i++) {
        for (int j = 0; j < _xSize; j++) {
            delete[] dist[i][j];
            delete[] degree[i][j];
        }
        delete[] dist[i];
        delete[] degree[i];
    }
    delete[] dist;
    delete[] degree;
}

void FCMOperator::initialization(int iNum, int cNum, int maxIter, double toler, double m) {
    imageNum = iNum;
    clusterNum = cNum;
    maxIteration = maxIter;
    tolerance = toler;
    weight = m;
}

void FCMOperator::inputLayer(vector<RasterLayer<double> *> layerD) {
    for (size_t i = 0; i < layerD.size(); ++i) {
        _vInputLayer.push_back(layerD[i]);
        Configure(layerD[i], false);
    }
    _noData = layerD[0]->_pMetaData->noData;
    _cellSize = layerD[0]->_pMetaData->cellSize;
    _xSize = layerD[0]->_pMetaData->_localdims.nRows();
    _ySize = layerD[0]->_pMetaData->_localdims.nCols();
    _nRows = layerD[0]->_pMetaData->row;
    _nCols = layerD[0]->_pMetaData->column;
    _rank = layerD[0]->_pMetaData->myrank;
    _pComptLayer = nullptr;
}

void FCMOperator::fcmLayer(RasterLayer<double>& layerD) {
    _pFCMLayer = &layerD;
    Configure(_pFCMLayer, false);
}

void FCMOperator::degLayer(vector<RasterLayer<double> *> layerD) {
    for (size_t i = 0; i < layerD.size(); ++i) {
        _vDegLayer.push_back(layerD[i]);
        Configure(layerD[i], false);
    }
}

bool FCMOperator::isTermination() {
    return flag;
}


bool FCMOperator::isNoData(double value, int rasterLayerIndex){
    return fabs(value - _vInputLayer[rasterLayerIndex]->metaData()->noData) <Eps;
}

bool FCMOperator::allNoDataAt(int row, int col) {
    for (int i=0;i<_vInputLayer.size();i++) {
        if(!isNoData((*_vInputLayer[i]->cellSpace())[row][col],i)) {
            return false;
        }
    }
    return true;
}
bool FCMOperator::anyNoDataAt(int row, int col) {
    for (int i=0;i<_vInputLayer.size();i++) {
        if(isNoData((*_vInputLayer[i]->cellSpace())[row][col],i)) {
            return true;
        }
    }
    return false;
}
inline bool FCMOperator::firstNoDataAt(int row, int col) {
    if(isNoData((*_vInputLayer[0]->cellSpace())[row][col],0)) {
        return true;
    }
    return false;
}
void FCMOperator::createRandomIdx(int nums, int range, int* randomIdx) {
    //generate nums numbers within the range
    srand((unsigned int)time(NULL)); //seed
    for (int i = 0; i < nums; i++) {
        int tmp = rand() % range; //0~n-1
        int j;
        for (j = 0; j < i; j++) {
            if (randomIdx[j] == tmp)
                break;
        }
        if (j >= i) {
            //nonredundant
            randomIdx[i] = tmp;
        }
        else {
            i--; //generate a new one if repeated
        }
    }
}

void FCMOperator::fnDistance(int curRow, int curCol, double* pInputVal) {
    for (int i = 0; i < clusterNum; i++) {
        dist[i][curRow][curCol] = 0.0;
        int count = 0;
        for (int j = 0; j < imageNum; j++) {
            //if (isNoData(pInputVal[j], j)) {
            //    cout<<"found a nodata! ("<< curRow<<","<<curCol<<")"<<endl;
            //    continue;
            //}
            dist[i][curRow][curCol] += (pInputVal[j] - centerVal[i * imageNum + j]) * (pInputVal[j] - centerVal[i *imageNum + j]);
            count++;//wyj 2020-10-21 update, calculate the mean value of the existing values
        }
        dist[i][curRow][curCol] = sqrt(dist[i][curRow][curCol]/count);
    }
}

//the membership degree from the cell to all clusters
void FCMOperator::InitDegree(int curRow, int curCol) {
    for (int p = 0; p < clusterNum; p++) {
        double sumDistance = 0.0;
        for (int q = 0; q < clusterNum; q++) {
            if (dist[q][curRow][curCol] == 0) {
                dist[q][curRow][curCol] = Eps;
            }
            sumDistance += pow(dist[p][curRow][curCol] / dist[q][curRow][curCol], 2 / (weight - 1));
        }
        degree[p][curRow][curCol] = sumDistance == 0.0 ? 1.0 : 1.0 / sumDistance;
        //first added to val for further judgement
        subval += pow(degree[p][curRow][curCol], weight) * pow(dist[p][curRow][curCol], 2);
    }
}

void FCMOperator::initRandomClusterCenters(double* clusterCenters) {
    cout << "initialization " << _noData << " is ok" << endl;
    //find all nonempty cells, and find clusterNum centers
    vector<int> tmpIdx;
    for (int i = 1; i < _xSize - 1; ++i) {
        for (int j = 1; j < _ySize - 1; ++j) {
            const double value = (*_vInputLayer[0]->cellSpace())[i][j];
            if (!isNoData(value,0)) {
                tmpIdx.push_back(i * _ySize + j);
            }
        }
    }
    int* randomIdx = new int[clusterNum];
    cout << "cells with value to init cluster center are " << tmpIdx.size() << endl;
    createRandomIdx(clusterNum, tmpIdx.size(), randomIdx); //put the center points into centerIndex
    cout << "init center idx and value is done. " << endl;
    for (int i = 0; i < clusterNum; ++i) {
        centerIndex[i] = tmpIdx[randomIdx[i]];
        int tmpRow = centerIndex[i] / _ySize;
        int tmpCol = centerIndex[i] % _ySize; 
        //get the center values
        for (int j = 0; j < imageNum; ++j) {
            centerVal[i * imageNum + j] = (*_vInputLayer[j]->cellSpace())[tmpRow][tmpCol];
        }
    }
}

///Clustering. find the maximum membership degree of every cell, and set the value to fcmL[i][j].
void FCMOperator::assignMaxMembershipDegrees() {
    CellSpace<double>& fcmL = *(_pFCMLayer->cellSpace());
    for (int i = 0; i < _xSize ; i++) {
        for (int j = 0; j < _ySize ; j++) {
            double time = MPI_Wtime();
            if (fabs((*_vInputLayer[0]->cellSpace())[i][j] - _noData) > Eps) {
            //if (!firstNoDataAt(i,j)) {
                int cNum = -1; //cluster index
                double degreeMax = 0.0;
                for (int p = 0; p < clusterNum; p++) {
                    if (degree[p][i][j] > degreeMax) {
                        degreeMax = degree[p][i][j];
                        cNum = p;
                    }
                    //to output degreeLayer
                    (*_vDegLayer[p]->cellSpace())[i][j] = degree[p][i][j];
                    //output entropy for test now.
                    if ((degree[p][i][j] - _vDegLayer[p]->metaData()->noData) > Eps) {
                        //this judgement is redundant
                        int nRows = _vInputLayer[0]->_pMetaData->row;
                        int nCols = _vInputLayer[0]->_pMetaData->column;
                        partitionCoef += degree[p][i][j] * degree[p][i][j] / (nRows * nCols - nodataNums);
                        if (degree[p][i][j] > 0) {
                            entropy += -1.0 * degree[p][i][j] * log(degree[p][i][j]) / (nRows * nCols -
                                nodataNums);
                        }
                    }
                }
                if(!_writePreExpLoad) {
                    fcmL[i][j] = cNum;
                }
            } else {
                if(!_writePreExpLoad){
                    fcmL[i][j] = _noData;
                }
                for (int p = 0; p < clusterNum; p++) {
                    (*_vDegLayer[p]->cellSpace())[i][j] = _noData;
                }
            }
            if (_writePreExpLoad){
                if (fcmL[i][j] < 0) {
                    fcmL[i][j] = 0;
                }
                fcmL[i][j] += (MPI_Wtime() - time) * 1000;
            } 
        }
    }
}

bool FCMOperator::Operator(const CellCoord& coord, bool operFlag) {
    int iRow = coord.iRow();
    int iCol = coord.iCol();
    CellSpace<double>& fcmL = *(_pFCMLayer->cellSpace());

    double* pInputVal = new double[imageNum];
    for (int i = 0; i < imageNum; ++i) {
        pInputVal[i] = (*_vInputLayer[i]->cellSpace())[iRow][iCol]; //[] first, then *
    }
    if ((iRow == 0) && (iCol == 0)) {
        if ((_iterNum == 0)) {
            //init for the first iteration
            centerVal = new double[clusterNum * imageNum]; //number * imageNum cluster centers
            centerIndex = new int[clusterNum]; // clusterNum * 1 
            sumNumerator = new double[clusterNum * imageNum];
            sumDenominator = new double[clusterNum];
            totNumerator = new double[clusterNum * imageNum];
            totDenominator = new double[clusterNum];
            //distance between points and cluster centers
            dist = new double **[clusterNum]; //dist 3-d, cluster numers * rows * cols
            for (int i = 0; i < clusterNum; ++i) {
                dist[i] = new double *[_xSize];
                for (int j = 0; j < _xSize; ++j) {
                    dist[i][j] = new double[_ySize];
                }
            }
            degree = new double **[clusterNum];
            for (int i = 0; i < clusterNum; i++) {
                degree[i] = new double *[_xSize];
                for (int j = 0; j < _xSize; j++) {
                    degree[i][j] = new double[_ySize];
                    for (int p = 0; p < _ySize; p++)
                        degree[i][j][p] = _noData;
                }
            }

            if (_rank == 0) {
                initRandomClusterCenters(centerVal);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(centerVal, clusterNum * imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //process 0 broadcast the cluster center
        }
        else {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(centerVal, clusterNum * imageNum, MPI_DOUBLE, 0, MPI_COMM_WORLD); //process 0 broadcast the cluster center
        }
    }
    
    //double time = MPI_Wtime();

    if (fabs(pInputVal[0] + 9999) <= Eps || fabs(pInputVal[0] - _noData) <= Eps) {
        
    }
    else {
        fnDistance(iRow, iCol, pInputVal); //distance from nonempty cells to the cluster centers
        InitDegree(iRow, iCol); //calculate membership degree
    }
    //FCMOperator::computeTimeExceptLastCell+=MPI_Wtime()-time;
    // if(allNoDataAt(iRow,iCol)) {
    //     for (int i = 0; i < clusterNum; i++) {
    //         dist[i][iRow][iCol] = -1;
    //     }
    // }else{
    //     fnDistance(iRow, iCol, pInputVal);
    //     InitDegree(iRow, iCol);
    // }

    //if (_writePreExpLoad){
    //    //(*_pComptLayer->cellSpace())[iRow][iCol] += (MPI_Wtime() - time) * 1000;//acutally not in use. Please use the "-out" file when setting "-dcmp compute -writeLoad /../../.."
    //    CellSpace<double>& fcmL = *(_pFCMLayer->cellSpace());
    //    if (fcmL[iRow][iCol] < 0) {
    //        fcmL[iRow][iCol] = 0;
    //    }
    //    fcmL[iRow][iCol] += (MPI_Wtime() - time) * 1000;
    //}
    
    //double time0 = MPI_Wtime();
    if ((iRow == _xSize - 1) && (iCol == _ySize - 1)) {
        //time = MPI_Wtime();
        MPI_Allreduce(&subval, &totval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //cout<<"rank"<<_rank<<" reduce1 time is "<<MPI_Wtime()-time<<"s"<<endl;
        _iterNum++;
        subval = 0.0;
        if (_rank == 0) {
            cout << "_iterNum " << _iterNum << " computation done. Start reducing." << endl;
        }
        //at last cell, val reduced to process 0, judge if the difference is under the threshold.
        if ((fabs(totval - oldtval) <= tolerance) || (_iterNum >= maxIteration)) {
            
            double tmpStart = MPI_Wtime();
            assignMaxMembershipDegrees();
            //time=MPI_Wtime();
            MPI_Allreduce(&partitionCoef, &totpartitionCoef, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&entropy, &totentropy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (_rank == 0) {
                double tmpEnd = MPI_Wtime();
                cout << "totpartitionCoef " << totpartitionCoef << " totentropy " << totentropy << " nodataNums " << nodataNums << endl;
                cout << _rank << " time is " << tmpEnd-tmpStart << endl;
            }
            //cout<<"rank"<<_rank<<" reduce2 time is "<<MPI_Wtime()-time<<"s"<<endl;
        }
        else {
            //time=MPI_Wtime();
            oldtval = totval;
            for (int p = 0; p < clusterNum; p++) {
                sumDenominator[p] = 0.0;
                for (int q = 0; q < imageNum; q++) {
                    sumNumerator[p * imageNum + q] = 0.0;
                }
                //find the new cluster center
                int valCount = 0;
                for (int i = 0; i < _xSize ; i++) {
                    for (int j = 0; j < _ySize ; j++) {
                        //double singleIterTime = MPI_Wtime();
                        //only calculate if the cell of first layer is nonempty.
                        if (fabs((*_vInputLayer[0]->cellSpace())[i][j] - _noData) > Eps && fabs((*_vInputLayer[0]->cellSpace())[i][j] + 9999) > Eps) {
                        //if (!firstNoDataAt(i,j)) {
                            sumDenominator[p] += pow(degree[p][i][j], weight);
                            for (int q = 0; q < imageNum; q++) {
                                sumNumerator[p * imageNum + q] += (pow(degree[p][i][j], weight) * (*_vInputLayer[q]->cellSpace())[i][j]);
                            }
                            ++valCount;
                        }
                        //if (_writePreExpLoad){
                        //    CellSpace<double>& fcmL = *(_pFCMLayer->cellSpace());
                        //    if (fcmL[i][j] < 0) {
                        //        fcmL[i][j] = 0;
                        //    }
                        //    fcmL[i][j] += (MPI_Wtime() - singleIterTime) * 1000;
                        //}
                        //multiple layer support
                        //if(allNoDataAt(i,j)) {
                        //    continue;
                        //}
                        //sumDenominator[p] += pow(degree[p][i][j], weight);
                        //for (int q = 0; q < imageNum; q++)
                        //{
                        //    const double value = (*_vInputLayer[q]->cellSpace())[i][j];
                        //    //if (isNoData(value,q)) {
                        //    //    continue;
                        //    //}
                        //    sumNumerator[p * imageNum + q] += pow(degree[p][i][j], weight) * value;
                        //    ++valCount;
                        //}
                    }
                }
            }
            //cout<<"rank"<<_rank<<" disp time is "<<MPI_Wtime()-time<<"s"<<endl;
            
            MPI_Allreduce(sumDenominator, totDenominator, clusterNum, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(sumNumerator, totNumerator, clusterNum * imageNum, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            //cout<<"rank"<<_rank<<" reduce3 time is "<<MPI_Wtime()-time<<"s"<<endl;
            for (int p = 0; p < clusterNum; p++) {
                for (int q = 0; q < imageNum; q++) {
                    centerVal[p * imageNum + q] = totNumerator[p * imageNum + q] / totDenominator[p];
                }
            }
            Termination = 0;
        }
    }
    //reduceTime+=MPI_Wtime()-time0;
    //FCMOperator::operatorReduceTime+=MPI_Wtime()-time0;

    delete pInputVal;

    return true;
}
