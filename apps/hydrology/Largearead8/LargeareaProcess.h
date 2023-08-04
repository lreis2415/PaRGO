#ifndef LARGEAREAPROCESS_H
#define LARGEAREAPROCESS_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include <cmath>
#include <functional>

using namespace std;
using namespace GPRO;

class LargeOperator : public RasterOperator<double> 
{
  public:
    LargeOperator()
      :RasterOperator<double>(),
       _pLayer(0),  _pwsLayer(0),_piLayer(0),_xSize(0), _ySize(0), num(0), flag(true),
	   srcfile(nullptr),outputfile(nullptr),nbr(nullptr),pfile(nullptr),wsfile(nullptr)
	{}
   
    ~LargeOperator() {}

	void getarea(RasterLayer<double> &layerD, int &minrow,int &mincol,int &nrow,int &ncol,int _g,int buf,int id);
	void init(RasterLayer<double>& layerD);
	void writeinto(RasterLayer<double> &layerD,RasterLayer<double> &layerD2,CellCoord nw);

	void initialfile();
	void writefile(CellCoord nw);
	CoordBR getsubarea(int &minrow,int &mincol,int &nrow,int &ncol, int _g,int buf,int id);

	//virtual bool isTermination();
    //virtual bool Operator(const CellCoord &coord,bool operFlag);

  public:
	double cellSize;
	double noData;
	int _xSize;	//rows in this processor
	int _ySize;
	int num;
	bool flag;
	int maxRow;
    int maxCol;
	int minRow;
	int minCol;
	char* outputfile;
	char* srcfile;
	char* wsfile;
	const char* pfile;
	char* nbr;
	RasterLayer<double> *_piLayer;
	RasterLayer<double> *_pwsLayer;
	RasterLayer<double> *_pLayer;
};
#endif