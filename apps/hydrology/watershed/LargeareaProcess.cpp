/***************************************************************************
* pPitRemove.cpp
*
* Project: GPRO_PitRemove
* Purpose: Pit&flat removing from DEM;Demonstration program for GPRO. 
*
* Author:  Ai Beibei
* E-mail:  aibb@lreis.ac.cn
****************************************************************************
* Copyright (c) 2015. Ai Beibei
* 
****************************************************************************/
#include"LargeareaProcess.h"

#define Eps 0.000001

void LargeOperator::init(RasterLayer<double>& layerD){	
	_piLayer=&layerD;
	CellSpace<double>& layer = *(_piLayer->cellSpace());
	
	maxRow = _piLayer->_pMetaData->_localworkBR.maxIRow();
    maxCol = _piLayer->_pMetaData->_localworkBR.maxICol();
	minRow = _piLayer->_pMetaData->_localworkBR.minIRow();
	minCol = _piLayer->_pMetaData->_localworkBR.minICol();
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			layer[i][j]=0;
		}	
	}
}
void LargeOperator::initialfile(){
	RasterLayer<double> Layer("Layer");
	Layer.readNeighborhood(nbr);  
	Layer.readFile(srcfile, NON_DCMP);
	RasterLayer<double> initLayer("initLayer");
	initLayer.copyLayerInfo(Layer);

	//LargeOperator LO;
	init(initLayer);
	initLayer.writeFile(outputfile);
}

void LargeOperator::getarea(RasterLayer<double>& layerD,int &minrow,int &mincol,int &nrow,int &ncol, int _g,int buf,int id) {
	_pwsLayer=&layerD;
	CellSpace<double>& watershed = *(_pwsLayer->cellSpace());

	int pi,pj,qi,qj;
	minRow = _pwsLayer->_pMetaData->_MBR.minIRow();
	minCol = _pwsLayer->_pMetaData->_MBR.minICol();
	maxRow = _pwsLayer->_pMetaData->_MBR.maxIRow();
	maxCol = _pwsLayer->_pMetaData->_MBR.maxICol();

	pi=-2;
	pj=-2;
	qi=-2;
	qj=-2;
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			if(watershed[i][j]==id){
				if(i<pi||pi==-2)
					pi=i;
				if(j<pj||pj==-2)
					pj=j;
				if(i>qi||qi==-2)
					qi=i;
				if(j>qj||qj==-2)
					qj=j;
			}
		}
	}
	
	if((pi-buf)>minRow)
		pi=pi-buf;
	else
		pi=minRow+1;
	if((pj-buf)>minCol)
		pj=pj-buf;
	else
		pj=minCol+1;
	if((qi+buf)<maxRow)
		qi=qi+buf;
	else
		qi=maxRow-1;
	if((qj+buf)<maxCol)
		qj=qj+buf;
	else
		qj=maxCol-1;

	minrow=pi*_g;
	mincol=pj*_g;
	nrow=(qi-pi+1)*_g;
	ncol=(qj-pj+1)*_g;
}
CoordBR LargeOperator::getsubarea(int &minrow,int &mincol,int &nrow,int &ncol, int _g,int buf,int id){
	RasterLayer<double> iLayer("iLayer");
	iLayer.readFile(outputfile,NON_DCMP);	
	RasterLayer<double> wsLayer("wsLayer");
	wsLayer.readFile(wsfile,NON_DCMP);

	getarea(wsLayer,minrow,mincol,nrow,ncol,_g,buf,id);
	
	//create block of sub-basin
	CellCoord nw(minrow,mincol);
	CellCoord se(minrow+nrow-1,mincol+ncol-1);
	CoordBR sub(nw,se);
	return sub;
}

void LargeOperator::writeinto(RasterLayer<double>& layerD,CellCoord nw){
	_pLayer=&layerD;
	CellSpace<double>& wk = *(_pLayer->cellSpace());
	CellSpace<double>& layer = *(_piLayer->cellSpace());
	cellSize = _pLayer->_pMetaData->cellSize;
    noData = _pLayer->_pMetaData->noData;
	//CellSpace<double>& dir = *(_pD8Layer->cellSpace());
	maxRow = _pLayer->_pMetaData->_localworkBR.maxIRow();
    maxCol = _pLayer->_pMetaData->_localworkBR.maxICol();
	int d1[9]={-1,-1,-1,0,0,0,1,1,1};
	int d2[9]={-1,0,1,-1,0,1,-1,0,1};
	int minRow = _pLayer->_pMetaData->_localworkBR.minIRow();
	int minCol = _pLayer->_pMetaData->_localworkBR.minICol();
	int gi = nw.iRow()-1;//_pwkLayer->_pMetaData->_MBR.minIRow();
	int gj = nw.iCol();//_pwkLayer->_pMetaData->_MBR.minICol();
	//cout<<"globali:"<<gi<<" globalj:"<<gj<<endl;
	//cout<<"scalayer"<<_pSCALayer->_pMetaData->_localworkBR.maxIRow()<<endl;
	int ti,tj;
	short d;
	double s=0.;
	for(int i=minRow;i<=maxRow;i++){
		for(int j=minCol;j<=maxCol;j++){
			//cout<<"sca:"<<sca[i+gi][j+gj]<<endl;
			//cout<<"wk:"<<wk[i][j]<<endl;
			if(fabs(layer[i+gi][j+gj]-noData)<=Eps||layer[i+gi][j+gj]<wk[i][j]){
				layer[i+gi][j+gj]=wk[i][j];
				//cout<<"wk:"<<wk[i][j]<<" sca:"<<sca[i+gi][j+gj]<<endl;
			}
			
				//cout<<"false"<<endl;

		}
	}
}
void LargeOperator::writefile(CellCoord nw){
	RasterLayer<double> pLayer("Layer");
	pLayer.readNeighborhood(nbr); 	 
	pLayer.readFile(pfile,NON_DCMP);
	RasterLayer<double> iLayer("iLayer");
	iLayer.readFile(outputfile,NON_DCMP);

	writeinto(pLayer,nw);
}