#include"LEOperator.h"
#include"math.h"
#include"geomorphons.h"
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include <iostream>
#include <fstream>

//Ϊ�˿���DP�ĵݹ����
/*directions
 * 3|2|1
 * 4|0|8
 * 5|6|7 */
//int nextr[8] = {-1, -1, -1, 0, 1, 1, 1, 0 };
//int nextc[8] = {1, 0, -1, -1, -1, 0, 1, 1 };

void LEOperator::demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  cellSize = _pDEMLayer->_pMetaData->cellSize;
  noData = _pDEMLayer->_pMetaData->noData;
    //pWorkBR = &_pDEMLayer->_pMetaData->_localworkBR;
  Configure(_pDEMLayer, false);
}

void LEOperator::LELayer(RasterLayer<double> &layerD) 
{
  _pLELayer = &layerD;
  Configure(_pLELayer,false);
}

bool LEOperator::isTermination()
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
    return true;
}

//�������
double caldist(int row1,int col1,int row2,int col2,double cell_res)
{
	double dist;
	dist = cell_res * sqrt((double)((col1-col2)*(col1-col2) +(row1-row2)*(row1-row2)));
	return dist;
}
//�����߲���
void ParaCalcu(vector<GeoPoint> Points,int firstpoint, int lastpoint,double *k,double *b)
{
	//LineParameters thisLine = new LineParameters();
	*k = ((double)Points[lastpoint].elevation-(double)Points[firstpoint].elevation)/((double)Points[lastpoint].distance-(double)Points[firstpoint].distance);
	*b = (double)Points[firstpoint].elevation - *k*(double)Points[firstpoint].distance;
	//cout<<"step1.1.1.1"<<" "<<firstpoint<<" "<<lastpoint<<endl;
	//cout<<"step1.1.1"<<" "<<Points[lastpoint].elevation<<" "<<Points[lastpoint].distance<<" "<<Points[firstpoint].elevation<<" "<<Points[firstpoint].distance<<endl;
	//cout<<"step1.1"<<" "<<k<<" "<<b<<endl;
}
//����ÿһ����ֱ�ߵĸ߲�
double DisCalcu(vector<GeoPoint> Points, int thispoint, double k,double b)
{
	//GeoPoint temPoint = new GeoPoint(0,0,0);
	double tem = 0; double distance;			
	tem = Points[thispoint].distance * k + b;			
	//distance = fabs(Points[thispoint].elevation - tem);
	double up = fabs((-k)*Points[thispoint].distance + Points[thispoint].elevation - b);
	//cout<<"up="<<up<<endl;
	double down = sqrt(k*k+1);
	//cout<<"down="<<down<<endl;
	distance = up/down;
	//cout<<"distance="<<distance<<endl;
	return distance;
}
//����DPһ�μ���
void DPCalcu(vector<GeoPoint> Points, int firstpoint, int lastpoint, double Tolerance, vector<int> &PointToKeep)
{
	double maxDistance = 0;
	int FarthestPoint = 0;
	double k,b;

	//����ÿһ����ֱ�ߵĸ߲�
	ParaCalcu(Points, firstpoint, lastpoint,&k,&b);
	//cout<<"step1"<<k<<" "<<b<<endl;
	for (int thispoint = firstpoint;thispoint<lastpoint;thispoint++)
	{
		double distance = DisCalcu(Points,thispoint,k,b);
		//cout<<"step2"<<thispoint<<" "<<distance<<endl;
		if (distance > maxDistance)
		{
			maxDistance = distance;
			FarthestPoint = thispoint;
		}
		//cout<<"step3"<<maxDistance<<" "<<FarthestPoint<<endl;
	}

		//��������㱣��
		PointToKeep.push_back(FarthestPoint);
		//�ݹ����
		//if (times == 0)
		//{
			//times = times + 1;
			//DPCalcu(Points, firstpoint, FarthestPoint, Tolerance,PointToKeep,times);
			//DPCalcu(Points, FarthestPoint, lastpoint, Tolerance, PointToKeep,times);
		//}
	//else
		return;
}

//����patten����
int determine_form(int num_minus, int num_plus)
{
/* determine form according number of positives and negatives
 * simple approach to determine form pattern */
	const FORMS forms[9][9] = 
	{
/* minus ------------- plus ----------------*/
/*       0   1   2   3   4   5   6   7   8  */
/* 0 */ {FL, FL, FL, FS, FS, VL, VL, VL, PT},
/* 1 */ {FL, FL, FS, FS, FS, VL, VL, VL, __},
/* 2 */ {FL, SH, SL, SL, CN, CN, VL, __, __},
/* 3 */ {SH, SH, SL, SL, SL, CN, __, __, __},
/* 4 */ {SH, SH, CV, SL, SL, __, __, __, __},
/* 5 */ {RI, RI, CV, CV, __, __, __, __, __},
/* 6 */ {RI, RI, RI, __, __, __, __, __, __},
/* 7 */ {RI, RI, __, __, __, __, __, __, __},
/* 8 */ {PK, __, __, __, __, __, __, __, __},
    };

/* legend:
  FL,  flat
  PK,  peak, summit
  RI,  ridge
  SH,  shoulder
  CV,  convex
  SL,  slope
  CN,  concave
  FS,  footslope
  VL,  valley
  PT,  pit, depression
  __  error, impossible
*/
 return forms[num_minus][num_plus];
}

int determine_form_num(int num_minus, int num_plus)
{
/* determine form according number of positives and negatives
 * simple approach to determine form pattern */
	const FORMS_NUM forms_num[9][9] = 
	{
		/* minus ------------- plus ----------------*/
		/*       0   1   2   3   4   5   6   7   8  */
		/* 0 */ {FLN, FLN, FLN, FSN, FSN, VLN, VLN, VLN, PTN},
		/* 1 */ {FLN, FLN, FSN, FSN, FSN, VLN, VLN, VLN, __N},
		/* 2 */ {FLN, SHN, SLN, SLN, CNN, CNN, VLN, __N, __N},
		/* 3 */ {SHN, SHN, SLN, SLN, SLN, CNN, __N, __N, __N},
		/* 4 */ {SHN, SHN, CVN, SLN, SLN, __N, __N, __N, __N},
		/* 5 */ {RIN, RIN, CVN, CVN, __N, __N, __N, __N, __N},
		/* 6 */ {RIN, RIN, RIN, __N, __N, __N, __N, __N, __N},
		/* 7 */ {RIN, RIN, __N, __N, __N, __N, __N, __N, __N},
		/* 8 */ {PKN, __N, __N, __N, __N, __N, __N, __N, __N},
	};

/* legend:
  FL,  flat
  PK,  peak, summit
  RI,  ridge
  SH,  shoulder
  CV,  convex
  SL,  slope
  CN,  concave
  FS,  footslope
  VL,  valley
  PT,  pit, depression
  __  error, impossible
*/
 return forms_num[num_minus][num_plus];
}

//�����ڼ���form
int extern multidetermin_form_num(int cur_Row, int cur_Col, vector<Pattern> patterns,int RowNum)
{
    vector<int> forms(RowNum-2);
	
	//�����д��ڵĽ��д��txt
	ofstream TXTout;
	TXTout.open("AllScalesLE.txt", ios::ate|ios::out|ios::app);
	if(TXTout.is_open())
	{
		TXTout<<cur_Row<<","<<cur_Col<<",";
	}
	
	for (int i=0;i<RowNum-2;i++)
	{
		forms[i] = determine_form_num(patterns[i].num_negatives,patterns[i].num_positives);	
		//����һ��formдһ��form
		if(TXTout.is_open())
		{
			TXTout<<forms[i]<<",";
		}		
	}
	//���դ��д���ˣ�Ҫд�����з�
	if(TXTout.is_open())
	{
		TXTout<<"\n";
	}

	//��¼����form���ִ���,���ϴ�������һ��11������
	int curform[11] = {0,0,0,0,0,0,0,0,0,0,0};
	//���������ڵ�form���м���
	
	for (int i = 0;i<RowNum-2;i++)
	{
		switch(forms[i])
		{
		case 1:curform[0]++;break;
		case 2:curform[1]++;break;
		case 3:curform[2]++;break;
		case 4:curform[3]++;break;
		case 5:curform[4]++;break;
		case 6:curform[5]++;break;
		case 7:curform[6]++;break;
		case 8:curform[7]++;break;
		case 9:curform[8]++;break;
		case 10:curform[9]++;break;
		default:curform[10]++;break;
		}
	}

	int formmax = 0;
	int formmaxnum = 0;
	for (int i=0;i<11;i++)
	{
		if (curform[i]>formmax)
		{
			formmax = i+1;
			formmaxnum = curform[i];
		}
	}
	int form = formmax;
	return form;
}
bool LEOperator::Operator(const CellCoord &coord, bool operFlag)
{
	//int rank;
	//MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//cout<<"rank"<<rank<<" ("<<coord.iRow()<<","<<coord.iCol()<<")"<<endl;
    CellSpace<double> &dem = *(_pDEMLayer->cellSpace());//����ͼ���դ������

    CellSpace<double> &LE = *(_pLELayer->cellSpace());//���ͼ���դ������

    Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);//���������ļ�

    int nextr[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
    int nextc[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };

    //coord�Ƿ������ڵ�����ԭ�㣬Ҳ������Ҫ�����¶ȵĵ�
    int cur_Row = coord.iRow();//��ǰ��row
    int cur_Col = coord.iCol();//��ǰ��cols

    double cell_res = 30;
    Pattern pattern;
    //Pattern patterns[RowNum-2];
    pattern.num_negatives = 0;
    pattern.num_positives = 0;
    pattern.negatives = 0;
    pattern.positives = 0;
    //���ýǶ���ֵ
    double flat_distance = cell_res * 3;//Ϊ�˱��ⴰ��������������е㽫��3��֮�⿪ʼ����
    double flat_threshold = 0.017;//�������һ����
    double flat_threshold_height = flat_distance*tan(flat_threshold);//�������һ����

    vector<GeoPoint> Points;//����㼯
    vector<int> PointToKeep;//���챣���㼯

    for (int i = 0; i<8; i++)
    {
        //�ര�ڼ��㣬�����Ҫ�ӵ�3������㵽��RowNum����
        int j = _icurrentScale;
        pattern.pattern[i] = 0;

        //��ʼ��Points
        GeoPoint curPoint;
        curPoint.i = i;
        for (int k = 0; k < j + 3; k++)
        {
            curPoint.row = cur_Row + k*nextr[i];
            curPoint.col = cur_Col + k*nextc[i];
            curPoint.elevation = dem[curPoint.row][curPoint.col];
            curPoint.height = curPoint.elevation - dem[cur_Row][cur_Col];
            curPoint.distance = cell_res * sqrt((double)((curPoint.col - cur_Col)*(curPoint.col - cur_Col) + (curPoint.row - cur_Row)*(curPoint.row - cur_Row)));
            curPoint.angel = atan2(curPoint.height, curPoint.distance);
            Points.push_back(curPoint);
        }

        //��ʼ�����������
        double Tolerance = 0;
        int firstpoint = 0;
        int lastpoint = Points.size() - 1;

        DPCalcu(Points, firstpoint, lastpoint, Tolerance, PointToKeep);

        int feature_point = PointToKeep[0];
        if (Points[feature_point].angel>flat_threshold)
        {
            pattern.positives += i;
            pattern.pattern[i] = 1;
            pattern.num_positives++;
            pattern.feature[i] = Points[feature_point];
        }
        if (Points[feature_point].angel<-flat_threshold)
        {
            pattern.negatives += i;
            pattern.pattern[i] = -1;
            pattern.num_negatives++;
            pattern.feature[i] = Points[feature_point];
        }
        else
        {
            pattern.pattern[i] = 0;
            pattern.feature[i] = curPoint;
        }
        //���ԭ��������
        Points.clear();
        PointToKeep.clear();
    }
    vector<GeoPoint>().swap(Points);
    vector<int>().swap(PointToKeep);
    int finalForm;
    finalForm = determine_form_num(pattern.num_negatives,pattern.num_positives);
    LE[cur_Row][cur_Col] = finalForm;
    return operFlag;
}