#include <ogrsf_frmts.h>
#include "utility.h"
#include "idwOperator.h"

inline double IDWOperator::getMinDistanceToBlockBound(double x,double y) {
    double boundXY[]={
        getBlockColIndexByCoord(x)*_blockSize+_sub_extent.minX, //left, block min x
        (getBlockColIndexByCoord(x)+1)*_blockSize+_sub_extent.minX, //right, block max x
        _sub_extent.maxY- getBlockRowIndexByCoord(y)*_blockSize, //down, block min y
        _sub_extent.maxY-(getBlockRowIndexByCoord(y)+1)*_blockSize //up, block max y
    };


    double distances[]={
        abs(x-boundXY[0]),
        abs(x-boundXY[1]),
        abs(y-boundXY[2]),
        abs(y-boundXY[3])
    };
    return min(min(distances[0],distances[1]),min(distances[2],distances[3]));
}

IDWOperator::~IDWOperator() {
    vector<SampleBlock> ().swap(_pSampleBlocks);
}

int IDWOperator::readSampleNums(const char* filename, char** pSpatialRefWkt) {
    //读取矢量样点的元数据，获取范围
#if GDAL_VERSION_MAJOR >= 2
	GDALAllRegister();
    GDALDataset* poDatasetsrc = (GDALDataset*)GDALOpenEx(filename, GDAL_OF_VECTOR, NULL, NULL, NULL);
#else
    OGRRegisterAll();
    OGRDataSource* poDatasetsrc = OGRSFDriverRegistrar::Open(filename, FALSE);
#endif
    if (poDatasetsrc == NULL) {
        printf("[ERROR] Open failed.\n");
        exit(1);
    }
    string file = filename;
    string f2 = file.substr(0, file.length() - 4);
    int pos = f2.find_last_of(SEP); //注意，linux用'/',windows用'\\'
    string f3 = f2.substr(pos + 1);

    OGRLayer* poLayer = poDatasetsrc->GetLayerByName(f3.c_str()); //f3是文件名，不带后缀
    OGRSpatialReference* sref = poLayer->GetSpatialRef();
    sref->exportToWkt(pSpatialRefWkt);

    OGRFeature* poFeature;

    poLayer->ResetReading();
    while ((poFeature = poLayer->GetNextFeature()) != NULL) {
        _sample_nums++;
        OGRFeature::DestroyFeature(poFeature);
    }
    //_sample_nums = poLayer->GetFeatureCount();	//为什么不直接用这个函数
#if GDAL_VERSION_MAJOR >= 2
	GDALClose(poDatasetsrc);
#else
    OGRDataSource::DestroyDataSource(poDatasetsrc);
#endif
    return _sample_nums;
}

bool IDWOperator::readSamples(const char* filename, int fieldIdx, char** pSpatialRefWkt, vector<SamplePoint> &samples) {

    //将位置信息和属性信息存放在数组Sample_Array中
#if GDAL_VERSION_MAJOR >= 2
	GDALAllRegister();
    GDALDataset* poDatasetsrc = (GDALDataset*)GDALOpenEx(filename, GDAL_OF_VECTOR, NULL, NULL, NULL);
#else
    OGRRegisterAll();
    OGRDataSource* poDatasetsrc = OGRSFDriverRegistrar::Open(filename, FALSE);
#endif
    if (poDatasetsrc == NULL) {
        printf("[ERROR] Open failed.\n");
        exit(1);
    }

    string file = filename;
    string f2 = file.substr(0, file.length() - 4);
    int pos = f2.find_last_of(SEP);
    string f3 = f2.substr(pos + 1);
    OGRLayer* poLayer = poDatasetsrc->GetLayerByName(f3.c_str());
    poLayer->ResetReading();

    int idx = 0;
    double x = 0.0;
    double y = 0.0;
    OGRFeature* poFeature;
    while ((poFeature = poLayer->GetNextFeature()) != NULL) {
        //cout<<"value:"<<showpoint<<Sample_Array[idx][2]<<endl;	//这里没问题，是小数
        OGRGeometry* poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint) {
            //读范围这一步为什么不用poLayer->GetExtent？
            OGRPoint* poPoint = (OGRPoint *)poGeometry;
            x = poPoint->getX();
            y = poPoint->getY();
            //存储位置信息
        	SamplePoint point;
        	point.x=x;
        	point.y=y;
        	point.value=poFeature->GetFieldAsDouble(fieldIdx);
        	samples.push_back(point);
            if (idx == 0) {
                _glb_extent.minX = x;
                _glb_extent.maxX = x;
                _glb_extent.minY = y;
                _glb_extent.maxY = y;
            }
            else {
                if (x > _glb_extent.maxX)
                    _glb_extent.maxX = x;
                if (x < _glb_extent.minX)
                    _glb_extent.minX = x;
                if (y > _glb_extent.maxY)
                    _glb_extent.maxY = y;
                if (y < _glb_extent.minY)
                    _glb_extent.minY = y;
            }
        }
        else {
            printf("[ERROR] No point geometry\n");
            return 1;
        }
        OGRFeature::DestroyFeature(poFeature);
        idx++;
    }

    _glb_extent.minX = _glb_extent.minX - _cellSize / 2;
    int totalcol = (_glb_extent.maxX - _glb_extent.minX) / _cellSize;
    if ((_glb_extent.maxX - _glb_extent.minX) != totalcol * _cellSize) {
        totalcol++;
    }
    _glb_extent.maxX = _glb_extent.minX + _cellSize * totalcol;

    _glb_extent.minY = _glb_extent.minY - _cellSize / 2;
    int totalrow = (_glb_extent.maxY - _glb_extent.minY) / _cellSize;
    if ((_glb_extent.maxY - _glb_extent.minY) != totalrow * _cellSize) {
        totalrow++;
    }
    _glb_extent.maxY = _glb_extent.minY + _cellSize * totalrow;
    //cout<<idw.extent_All.minX<<"	"<<idw.extent_All.maxX<<endl;
    _nRows = totalrow;
    _nCols = totalcol;

#if GDAL_VERSION_MAJOR >= 2
	GDALClose(poDatasetsrc);
#else
    OGRDataSource::DestroyDataSource(poDatasetsrc);
#endif
    return true;
}


void IDWOperator::creatSampleBlocks(vector<SamplePoint> &samples) {
    //思路：对idwLayer按行逐个分块，1D存储,获取该块的范围,暂时不需要范围
    //注意样点坐标范围和栅格有偏移；先求总行列号；
    //每个已知样点的xy都可推出所在块;(x-minX)/cellSize/blockGrain就可以求出行号,同理求列号；
    _blockRows = ceil((_glb_extent.maxY - _glb_extent.minY) / _blockSize);
    _blockCols = ceil((_glb_extent.maxX - _glb_extent.minX) / _blockSize);
    _pSampleBlocks.resize(_blockRows * _blockCols);
    for (int i = 0; i < _sample_nums; ++i) { 
        double x = samples[i].x;
        double y = samples[i].y; 
        int iCol = getBlockColIndexByCoord(x);
        int iRow = getBlockRowIndexByCoord(y);
        int index = iRow * _blockCols + iCol;
        if( index>= _pSampleBlocks.size()) {
            cerr<<"Index out of range. An error in creatSampleBlocks()"<<endl;
            cout<<"__blockRows * _blockCols = "<<_blockRows<<" * "<<_blockCols<<endl;
            printf("_pSampleBlocks.size()=%llu,samples.size()=%llu, _sample_nums=%i\n",_pSampleBlocks.size(),samples.size(),_sample_nums);
        }
        _pSampleBlocks[index].samplePoints.push_back(samples[i]);
    }
}
void IDWOperator::idwLayer(RasterLayer<double>& layerD, char** pSpatialRefWkt, DomDcmpType dcmpType) {
    //更新_pIDWLayer/layerD的基本元数据;即根据extent和_cellSize信息，创建栅格图层
    //MetaData **pMetaData = &(layerD._pMetaData);	//可以考虑用指针的指针简写
    layerD._pMetaData = new MetaData();
    if (layerD._pMetaData == NULL) {
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit(1);
    }

    layerD._pMetaData->noData = _noData;
    layerD._pMetaData->row = _nRows;
    layerD._pMetaData->column = _nCols;
    SpaceDims sdim(layerD._pMetaData->row, layerD._pMetaData->column);
    layerD._pMetaData->_glbDims = sdim;
    layerD._pMetaData->cellSize = _cellSize;
    layerD._pMetaData->format = "GTiff";
    layerD._pMetaData->_domDcmpType = ROWWISE_DCMP; //是否需要在这里指定
    MPI_Comm_rank(MPI_COMM_WORLD, &layerD._pMetaData->myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &layerD._pMetaData->processor_number);
    _myRank = layerD._pMetaData->myrank;

    DeComposition<double> deComp(layerD._pMetaData->_glbDims, *(layerD.nbrhood()));

    deComp.rowDcmp(*(layerD._pMetaData), layerD._pMetaData->processor_number); //根据数据范围按行划分,引用方式返回给*(layerD._pMetaData)

    layerD.newCellSpace(layerD._pMetaData->_localdims); //每个进程读入全区样点数据，但只算自己的workBR
    layerD._pMetaData->_glbDims.nRows(_nRows);
    layerD._pMetaData->_glbDims.nCols(_nCols);

    layerD._pMetaData->dataType = layerD.getGDALType();
    //pSpatialRefWkt 目前指向main函数中的char* pSpatialRefWkt的地址
    layerD._pMetaData->projection = *pSpatialRefWkt; //char* to string,直接赋值即可；string to char*,调用c_str()
    layerD._pMetaData->pTransform[0] = _glb_extent.minX;
    layerD._pMetaData->pTransform[1] = _cellSize;
    layerD._pMetaData->pTransform[2] = 0;
    layerD._pMetaData->pTransform[3] = _glb_extent.maxY;
    layerD._pMetaData->pTransform[4] = 0;
    layerD._pMetaData->pTransform[5] = -_cellSize;

    //更新子空间数据范围
    _sub_extent.minX = _glb_extent.minX;
    _sub_extent.maxX = _glb_extent.maxX;
    _sub_extent.maxY = _glb_extent.maxY - layerD._pMetaData->_MBR.minIRow() * _cellSize;
    _sub_extent.minY = _glb_extent.maxY - layerD._pMetaData->_MBR.maxIRow() * _cellSize - _cellSize;
    _xSize = layerD._pMetaData->_localdims.nCols();
    _ySize = layerD._pMetaData->_localdims.nRows();

    _pIDWLayer = &layerD;
    Configure(_pIDWLayer, false);
}
void IDWOperator::maskLayer(RasterLayer<int>& layerD) {
    _pMaskLayer = &layerD;
}
void IDWOperator::initIdwLayerGlobalInfo(RasterLayer<double>& layerD, char** pSpatialRefWkt) {
    //更新_pIDWLayer/layerD的基本元数据;即根据extent和_cellSize信息，创建栅格图层
    //MetaData **pMetaData = &(layerD._pMetaData);	//可以考虑用指针的指针简写
    layerD._pMetaData = new MetaData();
    if (layerD._pMetaData == NULL) {
        //do something
        cout << "[ERROR] MetaData is not allocate correct" << endl;
        exit(1);
    }

    layerD._pMetaData->noData = _noData;
    layerD._pMetaData->row = _nRows;
    layerD._pMetaData->column = _nCols;
    SpaceDims sdim(layerD._pMetaData->row, layerD._pMetaData->column);
    layerD._pMetaData->_glbDims = sdim;
    layerD._pMetaData->cellSize = _cellSize;
    layerD._pMetaData->format = "GTiff";
    layerD._pMetaData->_domDcmpType = ROWWISE_DCMP; //是否需要在这里指定
    MPI_Comm_rank(MPI_COMM_WORLD, &layerD._pMetaData->myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &layerD._pMetaData->processor_number);
    _myRank = layerD._pMetaData->myrank;

    layerD._pMetaData->_glbDims.nRows(_nRows);
    layerD._pMetaData->_glbDims.nCols(_nCols);

    layerD._pMetaData->dataType = layerD.getGDALType();
    //pSpatialRefWkt 目前指向main函数中的char* pSpatialRefWkt的地址
    layerD._pMetaData->projection = *pSpatialRefWkt; //char* to string,直接赋值即可；string to char*,调用c_str()
    layerD._pMetaData->pTransform[0] = _glb_extent.minX;
    layerD._pMetaData->pTransform[1] = _cellSize;
    layerD._pMetaData->pTransform[2] = 0;
    layerD._pMetaData->pTransform[3] = _glb_extent.maxY;
    layerD._pMetaData->pTransform[4] = 0;
    layerD._pMetaData->pTransform[5] = -_cellSize;

}
void IDWOperator::idwLayer(RasterLayer<double>& layerD, char** pSpatialRefWkt,CoordBR &subWorkBR) {

    initIdwLayerGlobalInfo(layerD,pSpatialRefWkt);


    layerD._pMetaData->_MBR=subWorkBR;
    layerD._pMetaData->_localdims=SpaceDims(subWorkBR.nRows(),subWorkBR.nCols());
    layerD.nbrhood()->calcWorkBR(layerD._pMetaData->_localworkBR,layerD._pMetaData->_localdims);
    layerD.newCellSpace(layerD._pMetaData->_localdims); //每个进程读入全区样点数据，但只算自己的workBR

    //更新子空间数据范围
    _sub_extent.minX = _glb_extent.minX;
    _sub_extent.maxX = _glb_extent.maxX;
    _sub_extent.maxY = _glb_extent.maxY - layerD._pMetaData->_MBR.minIRow() * _cellSize;
    _sub_extent.minY = _glb_extent.maxY - layerD._pMetaData->_MBR.maxIRow() * _cellSize - _cellSize;
    _xSize = layerD._pMetaData->_localdims.nCols();
    _ySize = layerD._pMetaData->_localdims.nRows();

    _pIDWLayer = &layerD;
    Configure(_pIDWLayer, false);
}
void IDWOperator::idwLayerSerial(RasterLayer<double>& layerD, char** pSpatialRefWkt) {
    if(GetRank()!=0) {
        return;
    }

    layerD._pMetaData = new MetaData();
    if (layerD._pMetaData == NULL) {
        cout << "[ERROR] MetaData not allocated." << endl;
        exit(1);
    }

    layerD._pMetaData->noData = _noData;
    layerD._pMetaData->row = _nRows;
    layerD._pMetaData->column = _nCols;
    SpaceDims sdim(layerD._pMetaData->row, layerD._pMetaData->column);
    layerD._pMetaData->_glbDims = sdim;
    layerD._pMetaData->cellSize = _cellSize;
    layerD._pMetaData->format = "GTiff";
    layerD._pMetaData->_domDcmpType = NON_DCMP;
    MPI_Comm_rank(MPI_COMM_WORLD, &layerD._pMetaData->myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &layerD._pMetaData->processor_number);
    _myRank = layerD._pMetaData->myrank;

    DeComposition<double> deComp(layerD._pMetaData->_glbDims, *(layerD.nbrhood()));

    deComp.rowDcmp(*(layerD._pMetaData), 1);

    layerD.newCellSpace(layerD._pMetaData->_localdims);
    layerD._pMetaData->_glbDims.nRows(_nRows);
    layerD._pMetaData->_glbDims.nCols(_nCols);

    layerD._pMetaData->dataType = layerD.getGDALType();
    layerD._pMetaData->projection = *pSpatialRefWkt;
    layerD._pMetaData->pTransform[0] = _glb_extent.minX;
    layerD._pMetaData->pTransform[1] = _cellSize;
    layerD._pMetaData->pTransform[2] = 0;
    layerD._pMetaData->pTransform[3] = _glb_extent.maxY;
    layerD._pMetaData->pTransform[4] = 0;
    layerD._pMetaData->pTransform[5] = -_cellSize;

    _sub_extent.minX = _glb_extent.minX;
    _sub_extent.maxX = _glb_extent.maxX;
    _sub_extent.maxY = _glb_extent.maxY - layerD._pMetaData->_MBR.minIRow() * _cellSize;
    _sub_extent.minY = _glb_extent.maxY - layerD._pMetaData->_MBR.maxIRow() * _cellSize - _cellSize;
    _xSize = layerD._pMetaData->_localdims.nCols();
    _ySize = layerD._pMetaData->_localdims.nRows();

    _pIDWLayer = &layerD;
    Configure(_pIDWLayer, false);
}
bool IDWOperator::isTermination() {
    return flag;
}

int IDWOperator::searchNbrSamples(const int subMinRow, int cellRow, int cellCol, double* nbrSamples) {
    double cellX = getXByCellIndex(cellCol); //当前待插值栅格坐标
    double cellY = getYByCellIndex(cellRow+subMinRow);
    int blockRow = getBlockRowIndexByCoord(cellY); //确定当前栅格所在块
    int blockRow1 = getBlockRowIndexByCellIndex(cellRow+subMinRow); //确定当前栅格所在块
    int blockCol = getBlockColIndexByCoord(cellX);
    int blockRows = _blockRows;
    int blockCols = _blockCols;

    double maxDist = 0.0; //目前搜索到的最大距离值;也许会受缓冲区限制
    int maxDistIdx = -1; //目前搜索到的最大距离样点所在位置
    int tailIdx = -1; //目前搜索到的样点尾部序列，即已搜索到的样点个数-1
    int searchRad = 0; //环形向外搜索半径;1代表3*3邻域
    bool isSearch = true;
    double minDistToBound=getMinDistanceToBlockBound(cellX,cellY);
    while (isSearch){
        double searchRange=(double)(searchRad-1) * _blockSize + minDistToBound;
        if(_idw_buffer>0 && searchRange >= _idw_buffer) {
            break;
        }

        //收集本层搜索的候选block idx
        vector<int> block2search;
        if (searchRad == 0) {
            block2search.resize(1);
        }
        else {
            block2search.resize(2 * searchRad * 4);
        }
        int blockCount = 0;
        for (int tRow = blockRow - searchRad; tRow <= blockRow + searchRad; ++tRow) {
            if (tRow < 0 || tRow >= blockRows) {
                continue;
            }
            if (tRow == blockRow - searchRad || tRow == blockRow + searchRad) {
                //首末两行存全部
                for (int tCol = blockCol - searchRad; tCol <= blockCol + searchRad; ++tCol) {
                    if (tCol < 0 || tCol >= blockCols) {
                        continue;
                    }
                    block2search[blockCount++] = tRow * blockCols + tCol;
                }
            }
            else {
                //其他行最多存左右边界两列,若边界两列有效则存，无效则跳过
                if (blockCol - searchRad >= 0 && blockCol - searchRad < blockCols) {
                    block2search[blockCount++] = tRow * blockCols + blockCol - searchRad;
                }
                if (blockCol + searchRad >= 0 && blockCol + searchRad < blockCols) {
                    block2search[blockCount++] = tRow * blockCols + blockCol + searchRad;
                }
            }
        }
        //cout<<"myrank "<<_myRank<<" "<<tRow<<" "<<tCol<<" "<<_pSampleBlocks[tRow*blockCols+tCol].samplePoints.size()<<" "<<endl;
        //遍历int block2search[2*searchRad*4]中存储的block idx;
        //对_pSampleBlocks[i].samplePoints进行搜索
        for (int i = 0; i < blockCount; ++i) {
            int blockIdx = block2search[i];
            //cout<<blockIdx<<" "<<_pSampleBlocks[blockIdx].samplePoints.size()<<endl;
            for (vector<SamplePoint>::iterator iter = _pSampleBlocks[blockIdx].samplePoints.begin(); iter != _pSampleBlocks[blockIdx].samplePoints.end(); ++iter) {
                double tmpDist = sqrt(pow(iter->x - cellX,2) + pow(iter->y - cellY,2));
                if(tmpDist==0) {
                    tmpDist=EPS;
                }

                if (tailIdx < _nbrPoints - 1) {
                    //搜索到的点还不足指定的个数，则直接放入尾部
                    tailIdx++;
                    nbrSamples[tailIdx * 2] = tmpDist;
                    nbrSamples[tailIdx * 2 + 1] = iter->value;
                    if (tmpDist > maxDist) {
                        maxDist = tmpDist;
                        maxDistIdx = tailIdx;
                    }
                }
                else {
                    tailIdx++;
                    //已经有足够邻近样点，若所搜到的新点距离更近，则替换目前最远那个点,并更新最远距离及ID
                    if (tmpDist < maxDist) {
                        nbrSamples[maxDistIdx * 2] = tmpDist;
                        nbrSamples[maxDistIdx * 2 + 1] = iter->value;
                        maxDist = nbrSamples[0];
                        maxDistIdx = 0;
                        //更新maxDist,考虑改用有序数据结构，即可省去这里
                        for (int i = 1; i < _nbrPoints; ++i) {
                            if (nbrSamples[i * 2] > maxDist) {
                                maxDist = nbrSamples[i * 2];
                                maxDistIdx = i;
                            }
                        }
                    }
                }
            }
        }
        if (tailIdx >= _nbrPoints - 1) {
            if (searchRange >= maxDist) {
                isSearch=false;
            }
            else {
                ++searchRad;
            }
        }
        else {
            ++searchRad;
        }
        //delete []block2search;
        //block2search=nullptr;
    }
    return min(tailIdx+1,_nbrPoints);
}

bool IDWOperator::Operator(const CellCoord& coord, bool operFlag) {
    double startTime = MPI_Wtime();
    int iRow = coord.iRow();
    int iCol = coord.iCol();

    int maskRow=_pIDWLayer->rowAtOtherLayer(_pMaskLayer,iRow);
    int maskCol=_pIDWLayer->colAtOtherLayer(_pMaskLayer,iCol);
    double mask = (*_pMaskLayer->cellSpace())[maskRow][maskCol];
    double maskNoData=_pMaskLayer->metaData()->noData;
    if(mask==maskNoData) {
        (*_pIDWLayer->cellSpace())[iRow][iCol]=_noData;
        if(_pComptLayer) {
            (*_pComptLayer->cellSpace())[iRow][iCol] += (MPI_Wtime()-startTime) * 1000;
        }
        return true;
    }

    int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    CellSpace<double>& idwL = *_pIDWLayer->cellSpace();
    const int minRow = _pIDWLayer->_pMetaData->_MBR.minIRow();
 
    //每个点都是待插值点，只是搜索范围不同而已
    double* pNbrSamples = new double [_nbrPoints * 2]; //依次存放距离和属性值对
    int sampleNum=searchNbrSamples(minRow, iRow, iCol, pNbrSamples); //搜索当前栅格的样点值，在nbrSamples中返回
    //
    ////计算插值结果
    double weightSum = 0.0;
    double* pWeight = new double[_nbrPoints];
    idwL[iRow][iCol]=0;
    for (int i = 0; i < sampleNum; ++i) {
        pWeight[i] = 1 / pow(pNbrSamples[i * 2], _idw_power);
        weightSum += pWeight[i];
    }
    for (int i = 0; i < sampleNum; ++i) {
        double value=pNbrSamples[i * 2 + 1] * pWeight[i] / weightSum;
        //idwL[iRow][iCol] += pNbrSamples[i * 2 + 1] * pWeight[i] / weightSum;
        idwL[iRow][iCol] += value;
    }
	
//for test.
    //idwL[iRow][iCol] = minRow;
    //int blockRow = getBlockRowIndexByCellIndex(iRow+minRow);
    //int blockCol = getBlockColIndexByCellIndex(iCol);
    //int blockCols = _blockCols;
    ////idwL[iRow][iCol] = blockRow*blockCols+blockCol;
    //idwL[iRow][iCol] = 100;

    if(_pComptLayer) {
        (*_pComptLayer->cellSpace())[iRow][iCol] += (MPI_Wtime()-startTime) * 1000;
    }
    delete []pNbrSamples;
    pNbrSamples=nullptr;
    delete []pWeight;
    pWeight=nullptr;
    
    return true;
}
