#include"d8Operator.h"
//#include "math.h"
#define Eps 0.000001

void D8Operator::
demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  cellSize = _pDEMLayer->_pMetaData->cellSize;
  noData = _pDEMLayer->_pMetaData->noData;
  Configure(_pDEMLayer, false);
}

void D8Operator::d8Layer(RasterLayer<double> &layerD) 
{
  _pD8Layer = &layerD;
  Configure(_pD8Layer,false);
}

bool D8Operator::isTermination()
{
	num--;
    return num > 0;
}


bool D8Operator::Operator(const CellCoord &coord,bool operFlag)
{
	CellSpace<double> &demL = *(_pDEMLayer->cellSpace());
	CellSpace<double> &d8L = *(_pD8Layer->cellSpace());
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);
	
	int iRow = coord.iRow();
	int iCol = coord.iCol();
	int iNeighborCells = ((int)sqrt((double)nbrhoodD.size())) / 2;
	double d[9];
	int k = 0;
	bool tag = true;

	for(int iRow1 = iRow - iNeighborCells; iRow1 <= iRow + iNeighborCells; iRow1++)
	{
		for(int iCol1 = iCol - iNeighborCells; iCol1 <= iCol + iNeighborCells; iCol1++)
		{
			d[k] = demL[iRow1][iCol1];
			if( fabs(d[k]-noData)<Eps )
			{
				tag = false;
			}
			k++;
		}
	}
	if( !tag ){
		d8L[iRow][iCol] = noData;
		return true;
	}
	else
	{
		//i��������max{k*(Zc-Zi)};i=1,2������8
		//��iλ�ڶ����ϱ�ʱ����k=1;��iΪ�Խ��߷���ʱ��k=1/1.414
		int valid[9] = {0};
		double max = 0.0;
		double diff[9];
		for(int i=0;i<9;i++)
		{
			diff[i] = d[4]-d[i];
			if( i==0||i==2||i==6||i==8 )
				diff[i] = diff[i]/sqrt(2.0);
			//      diff[i]= diff[i]/1.5;
			if(diff[i]>max)
				max = diff[i];
		}
		if( max < Eps ){
			d8L[iRow][iCol] = -1;
			return true;
		}
		//����dirType��֧ѡ�㷨
		if( dirType == 8 ){
			//���ֻ��0-7����ͬ���½�����������ѡ�м��һ����ż����ѡ�м俿���һ��
			vector<int> validDir;
			for(int i=0;i<9;i++)
			{
				if( fabs(diff[i] - max)<Eps ){
					validDir.push_back(i);
				}
			}
			int length = validDir.size();
			d8L[iRow][iCol] = validDir[length>>1];
		}else{
			//��ͬ���½������ۼӷ���ֵ�����õ�������뷽ʽͬarcgis10.0
			for(int i=0;i<9;i++)
			{
				if( fabs(diff[i] - max)<Eps )
					valid[i] = 1;
			}
			d8L[iRow][iCol] = valid[5]*pow(2.0,0)+valid[8]*pow(2.0,1)+valid[7]*pow(2.0,2)+valid[6]*pow(2.0,3)
								+valid[3]*pow(2.0,4)+valid[0]*pow(2.0,5)+valid[1]*pow(2.0,6)+valid[2]*pow(2.0,7);
		}
		return true;
	}
}

