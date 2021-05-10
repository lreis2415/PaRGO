#ifndef KRIGINGOPERATOR_H
#define KRIGINGOPERATOR_H

#include "cellSpace.h"
#include "neighborhood.h"
#include "rasterOperator.h"
#include "rasterLayer.h"
#include "matrix.h"
#include <cmath>
#include <functional>

#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>

#include <set>
#include <string>
#include <sstream>
#include <utility>
#include <vector>


using namespace GPRO;

#define EPS 0.0000001
#define  NODATA_DEFINE -9999




struct SamplePoint {
	double x;
	double y;
	double value;
	SamplePoint() {
	}
	SamplePoint(const SamplePoint& sourcePoint) {
		x = sourcePoint.x;
		y = sourcePoint.y;
		value = sourcePoint.value;
	}
};

struct extent_info {
	double minX;
	double maxX;
	double minY;
	double maxY;
};

struct SampleBlock {
	vector<SamplePoint> samplePoints;
};


class KrigingOperator : public RasterOperator<double>
{
public:
	KrigingOperator()
		:RasterOperator<double>(),
		_iterNum(0), flag(true), _sample_nums(0), _pMaskLayer(nullptr), estimatedNugget(0),estimatedSill(0),estimatedRange(0), lag(0),lagTolerance(0)
	{}

	KrigingOperator(float cellsize, int nbrPoints, double bufferSize, double blockSize)
		:RasterOperator<double>(),
		_iterNum(0), flag(true), _sample_nums(0), _pMaskLayer(nullptr), estimatedNugget(0),estimatedSill(0),estimatedRange(0), lag(0),lagTolerance(0), _cellSize(cellsize), _nbrPoints(nbrPoints), _idw_buffer(bufferSize), _blockSize(blockSize), _noData(NODATA_DEFINE)
	{}

	~KrigingOperator();

	double getBlockSize() { return _blockSize; }

	int readSampleNums(const char* filename, char** pSpatialRefWkt);
	bool readSamples(const char* filename, int fieldIdx, char** pSpatialRefWkt, vector<SamplePoint> &samples);
	void creatSampleBlocks(vector<SamplePoint> &samples);
	const vector<SampleBlock>* getSampleBlocks() { return &_pSampleBlocks; }

	void idwLayer(RasterLayer<double> &layerD, char** pSpatialRefWkt, DomDcmpType dcmpType = ROWWISE_DCMP);
	void idwLayer(RasterLayer<double> &layerD, char** pSpatialRefWkt, CoordBR& subWorkBR);
	void maskLayer(RasterLayer<int> &layerD);
	void maskLayer(RasterLayer<int> &layerD, CoordBR& subWorkBR);
	void idwLayerSerial(RasterLayer<double> &layerD, char** pSpatialRefWkt);
	virtual bool isTermination();

	std::multimap<double, SamplePoint> searchNbrSamples(const int subMinRow, int cellRow, int cellCol, double cellX, double cellY);

	virtual bool Operator(const CellCoord &coord, bool operFlag);

	int getBlockRowIndexByCoord(double y) { return (_glb_extent.maxY - y) / _blockSize; }
	int getBlockColIndexByCoord(double x) { return (x - _glb_extent.minX) / _blockSize; }
	int getBlockRowIndexByCellIndex(int iRow, double cellSize) { return min(int((iRow + 0.5)*cellSize / _blockSize), _blockRows - 1); }
	int getBlockColIndexByCellIndex(int iCol, double cellSize) { return min(int((iCol + 0.5)*cellSize / _blockSize), _blockCols - 1); }
	double getYByCellIndex(int iRow, double cellSize) { return double(_glb_extent.maxY - (iRow + 0.5) * cellSize); }
	double getXByCellIndex(int iCol, double cellSize) { return double((iCol + 0.5) * cellSize + _glb_extent.minX); }

	int getBlockCols() { return _blockCols; }
	int getBlockRows() { return _blockRows; }
	float getCellSize() { return _cellSize; }
	double getGlobalMinX() { return _glb_extent.minX; }
	double getGlobalMaxX() { return _glb_extent.maxX; }
	double getGlobalMinY() { return _glb_extent.minY; }
	double getGlobalMaxY() { return _glb_extent.maxY; }
	RasterLayer<int>* getMaskLayer() { return _pMaskLayer; }

	inline double getMinDistanceToBlockBound(double x, double y);
	int getNbrPoints() { return _nbrPoints; }
	void initIdwLayerGlobalInfo(RasterLayer<double>& layerD, char** pSpatialRefWkt);



	/******************
	* KRIGING MODULE
	*******************/
	enum Model
	{
		Linear,
		LinearWithoutIntercept,
		Spherical,
		Exponential,
		Gaussian
	};

	void Initialize();
	void ParallelInitialize();
	void CalculateExperimentalVariogram(double lag, double lagTolerance);

	std::vector<double> CalculateLinearModel() const;
	std::vector<double> CalculateLinearModelWithoutIntercept(double nugget) const;
	std::vector<double> CalculateSphericalModel(double nugget, double sill, double range) const;
	std::vector<double> CalculateExponentialModel(double nugget, double sill, double range) const;
	std::vector<double> CalculateGaussianModel(double nugget, double sill, double range) const;


	// Retrieve semivariogram parameters

	double GetEstimatedNugget() const
	{
		return estimatedNugget;
	}

	double GetEstimatedRange() const
	{
		return estimatedRange;
	}

	double GetEstimatedSill() const
	{
		return estimatedSill;
	}

	// Get lag parameters

	double GetDefaultLag() const
	{
		return lag;
	}

	double GetDefaultLagTolerance() const
	{
		return lagTolerance;
	}

	// Get empirical semivariogram vectors

	const std::vector<double>& GetLagDistances()
	{
		return lagDistance;
	}

	const std::vector<double>& GetLagSemivariances()
	{
		return lagSemivariance;
	}

	const std::vector<unsigned int>& GetLagCounts()
	{
		return lagCount;
	}

	// Ordinary Kriging

	void OrdinaryKrige(Model model, double nugget, double sill, double range,
		unsigned int minPoints, unsigned int maxPoints, double maxDistance);

private:
	int getBlockRowIndexByCellIndex(int iRow) { return (iRow + 0.5)*_cellSize / _blockSize; }
	int getBlockColIndexByCellIndex(int iCol) { return (iCol + 0.5)*_cellSize / _blockSize; }

	double getYByCellIndex(int iRow) { return double(_glb_extent.maxY - (iRow + 0.5) * _cellSize); }
	double getXByCellIndex(int iCol) { return double((iCol + 0.5) * _cellSize + _glb_extent.minX); }
    
    vector<int> decomposePointsForVariogram(int processNums, vector<int>& pairsPerProcess);
    int SFi(int i0, int i1){return ((_sample_nums-i0-1)+(_sample_nums-i1-1))*(i1-i0+1)/2;}
	float _cellSize;
	int _xSize, _ySize;
	int _nRows, _nCols;
	double _noData;
	int _myRank;
	int _iterNum;
	bool flag;
protected:
	int _nbrPoints;
	int _idw_power;
	double _idw_buffer;
	extent_info _glb_extent;
	extent_info _sub_extent;
	int _sample_nums;
	double _blockSize;
	vector<SampleBlock> _pSampleBlocks;
	int _blockRows;
	int _blockCols;
	RasterLayer<double> *_pIDWLayer;
	RasterLayer<int> *_pMaskLayer;

	void CalculatePartialDistanceMap(int start, int end, double* distances, double* variances);
	void CalculateDistanceMap();

	std::vector<std::vector<double>> CalculateVariogramMatrix(const std::vector<SamplePoint>& dataPointCandidate,Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;
	std::vector<std::vector<double>> CalculateCovariogramMatrix(const std::vector<SamplePoint>& dataPointCandidate,Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;
	std::vector<double> CalculateVariogramVector(const std::vector<SamplePoint>& dataPointCandidate,double xCoordinate, double yCoordinate, Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;
	std::vector<double> CalculateCovariogramVector(const std::vector<SamplePoint>& dataPointCandidate,double xCoordinate, double yCoordinate, Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const;
	std::multimap<double, size_t> CalculateDistanceMapForPoint(double pointX, double pointY, unsigned int maxPoints, double maxDistance) const;

	void CalculateDefaultLagParameters();

	void CalculateEstimatedVariogramParameters();

	double CalculateVariogram(Model model, double distance, double nugget, double sill, double range) const;

	double CalculateCovariogram(Model model, double distance, double nugget, double sill, double range) const;

	void DoSimpleLinearRegression(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double* intercept) const;
	void DoSimpleLinearRegressionWithoutIntercept(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double intercept) const;

	// Ordinary Kringing for an individual point

	double OrdinaryKrigeForPoint(double xCoordinate, double yCoordinate,
		Model model, double nugget, double sill, double range,
		LUDecomposition& luDecomposition, const std::vector<SamplePoint>& dataPointCandidate);

	double estimatedNugget;
	double estimatedSill;
	double estimatedRange;

	double lag;
	double lagTolerance;

	vector<SamplePoint> dataPoint;

	multimap<double, double> semiVariogram;

	vector<double> lagDistance;
	vector<double> lagSemivariance;
	vector<unsigned int> lagCount;
};


#endif