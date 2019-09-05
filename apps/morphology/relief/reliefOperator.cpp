#include"reliefOperator.h"


//添加输入图层
void reliefOperator::
demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  cellSize = _pDEMLayer->_pMetaData->cellSize;
  noData = _pDEMLayer->_pMetaData->noData;
    //pWorkBR = &_pDEMLayer->_pMetaData->_localworkBR;
  Configure(_pDEMLayer, false);
}

//添加输出图层
void reliefOperator::reliefLayer(RasterLayer<double> &layerD) 
{
  _pReliefLayer = &layerD;
  Configure(_pReliefLayer,false);
}

//判断迭代
bool reliefOperator::isTermination()
{
	num--;
	if(num > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//计算起伏度的算法实现
bool reliefOperator::Operator(const CellCoord &coord)
{
	CellSpace<double> &dem = *(_pDEMLayer->cellSpace());//输入图层的栅格数据

	CellSpace<double> &relief = *(_pReliefLayer->cellSpace());//输出图层的栅格数据
	  
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);//分析窗口
	
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	
	double d[9];//存放3*3邻域窗口
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	int dCellSize = _pDEMLayer->_pMetaData->cellSize;//DEM格网长度
	int nodata = _pDEMLayer->_pMetaData->noData;
	int iRow1, iCol1;
	
	//存储分析窗口范围内的DEM值
	int k = 0;
	int tag=0;
	for(iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++)
	{
		for(iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++)
		{
			d[k] = dem[iRow1][iCol1];
			if(d[k]==nodata)
			{
				tag=1;
			}
			k++;
		}
	}
			if(tag==1)
		{
			relief[iRow][iCol] = nodata;
			return true;
		}
		else
		{
			//起伏度（高差）算法
			double max,min;
			max=d[0];	min=d[0];
			for(int m=1;m<9;m++)
			{
				if(d[m]>max)
					max=d[m];
				else if(d[m]<min)
					min=d[m];
			}
			relief[iRow][iCol] =max-min;
			return true;
		}
}

