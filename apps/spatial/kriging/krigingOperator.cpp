#include <ogrsf_frmts.h>
#include "utility.h"
#include "krigingOperator.h"


void KrigingOperator::Initialize()
{
	CalculateDistanceMap();
	CalculateDefaultLagParameters();
	CalculateExperimentalVariogram(GetDefaultLag(), GetDefaultLagTolerance());
	CalculateEstimatedVariogramParameters();
    
    //cout<<"nugget: "<<GetEstimatedNugget()<<endl;
    //cout<<"sill: "<<GetEstimatedSill()<<endl;
    //cout<<"range: "<<GetEstimatedRange()<<endl;
    //cout<<"lag distance: "<<GetLagDistances()<<endl;
}

// Not in use. The reducing contains some parallelism details, and not much improvement is made
//   unless the sorting of distance-variance pairs is parallelized, which is another difficulty.
void KrigingOperator::ParallelInitialize()
{

	int processNums;
	MPI_Comm_size(MPI_COMM_WORLD, &processNums);
    vector<int> pairsEveryProcess;
    vector<int> pointPartitions=decomposePointsForVariogram(processNums,pairsEveryProcess);
    if (pointPartitions.size()<processNums) {
        cerr<< "decomposePointsForVariogram() failed to decompose the points into specified partitions"<<endl;
        exit(0);
    }
	int startIndex,endIndex;
	if (GetRank() == 0) {
		startIndex = 0;
	}
	else {
		startIndex = pointPartitions[GetRank()-1]+1;
	}
    endIndex = pointPartitions[GetRank()];
    int pointsInThisProcess = endIndex - startIndex+1;
	int pairsTotal = _sample_nums *(_sample_nums - 1) / 2;

	int pairsInThisProcess = SFi(startIndex,endIndex);

	double *distances = new double[pairsInThisProcess];
	double *variences = new double[pairsInThisProcess];

    
    double startTime = MPI_Wtime();
	CalculatePartialDistanceMap(startIndex,endIndex,distances,variences);
    cout<<"CalculatePartialDistanceMap done. time:"<<MPI_Wtime()-startTime<<endl;

    vector<int> displs;
    displs.emplace_back(0);
    for (int i = 0; i < pairsEveryProcess.size()-1; ++i) {
        displs.emplace_back(displs.back()+pairsEveryProcess[i]);
    }
	double *distancesRecv = new double[pairsTotal];
	double *variencesRecv = new double[pairsTotal];
	double sill=0;
    startTime = MPI_Wtime();
	MPI_Allgatherv(distances, pairsInThisProcess, MPI_DOUBLE, distancesRecv, &pairsEveryProcess[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgatherv(variences, pairsInThisProcess, MPI_DOUBLE, variencesRecv, &pairsEveryProcess[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allreduce(&estimatedSill, &sill, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    cout<<"Allgatherv time:"<<MPI_Wtime()-startTime<<endl;
	estimatedSill = (sill / 2.0) / ((dataPoint.size() * (dataPoint.size() - 1)) / 2.0);
    
    //for (int i = 0; i < pairsTotal; ++i) {
    //    if(distances[i]!=distancesRecv[i] || distancesRecv[i] <0 || distances[i]<0)
    //        cout<<i<<": "<<distances[i]<<" - "<<distancesRecv[i]<<endl;
    //}
    
	delete[]distances;
	distances = nullptr;
	delete[]variences;
	variences = nullptr;

    startTime = MPI_Wtime();
    for (int i = 0; i < pairsTotal; ++i) {
		semiVariogram.insert(pair<double, double>(distancesRecv[i], variencesRecv[i]));
    }
    cout<<"construct semivariogram of "<<pairsTotal<<" pairs multimap time:"<<MPI_Wtime()-startTime<<endl;
    
	delete[]distancesRecv;
	distancesRecv = nullptr;
	delete[]variencesRecv;
	variencesRecv = nullptr;
    
    cout<<"semiVariogram done."<<endl;
    
	CalculateDefaultLagParameters();

	CalculateExperimentalVariogram(GetDefaultLag(), GetDefaultLagTolerance());

	CalculateEstimatedVariogramParameters();
    
    startTime = MPI_Wtime();
    semiVariogram.clear();
    cout<<"clear semiVariogram multimap time."<<MPI_Wtime()-startTime<<endl;

    if (GetRank()==0) {
        cout<<"nugget: "<<GetEstimatedNugget()<<endl;
        cout<<"sill: "<<GetEstimatedSill()<<endl;
        cout<<"range: "<<GetEstimatedRange()<<endl;
        cout<<"lag distance: "<<GetLagDistances()<<endl;
    }
}

//return partitions, out pairsPerProcess
vector<int> KrigingOperator::decomposePointsForVariogram(int processNums, vector<int>& pairsPerProcess) {
    int rank=GetRank();
    vector<int> partitions;
    int k=processNums;
    int a=_sample_nums;
    double mean = double(a*(a-1))/(2*k);
    int p=0; // partition
    int sfiP=0;
    for (int i = 0; i < processNums-1; ++i) {
        long temp1t1 = a * a * k;
        long temp1t2 = a * a + 2 * a * p * k + 2 * a * k;
        long temp1t3 = a + p * p * k + 2 * p * k + k;
        long temp1sum = temp1t1 - temp1t2 + temp1t3;
        double temp1=sqrt(static_cast<double>(k) * temp1sum);
        //double temp1=sqrt(static_cast<double>(k) * (a * a * k - a * a - 2 * a * p * k - 2 * a * k + a + p * p * k + 2 * p * k + k));
        double temp2=a*k-p*k-k;
        int p1 = (-temp1 + temp2)/k;
        int p2 = (temp1 + temp2)/k;
        int deltaP = p1;
        if(p+p1<=0 || p+p1 >=_sample_nums) {
            deltaP = p2;
        }
        int sfiP1 = SFi(p,p+deltaP);
        int sfiP2 = SFi(p,p+deltaP+1);
        if(abs(sfiP1-mean)<abs(sfiP2-mean)) {
            p+=deltaP;
            sfiP=sfiP1;
        }else {
            p+=deltaP+1;
            sfiP=sfiP2;
        }
        partitions.emplace_back(p);
        pairsPerProcess.emplace_back(sfiP);
        p++;
    }
    partitions.emplace_back(_sample_nums-1);
    pairsPerProcess.emplace_back(SFi(p,_sample_nums-1));
    return partitions;
}

// Fill a map of distances from a fixed point with variogram

std::multimap<double, size_t> KrigingOperator::CalculateDistanceMapForPoint(double pointX, double pointY, unsigned int maxPoints, double maxDistance) const
{
	std::multimap<double, size_t> pointDistanceMap;

	double currentDistance = 0.0;

	for (size_t i = 0; i < dataPoint.size(); i++)
	{
		currentDistance = std::sqrt(std::pow(dataPoint[i].x - pointX, 2.0) + std::pow(dataPoint[i].y - pointY, 2.0));

		if (currentDistance <= maxDistance)
		{
			if (pointDistanceMap.size() == maxPoints)
			{
				if (currentDistance < pointDistanceMap.rbegin()->first)
				{
					pointDistanceMap.insert(std::pair<double, size_t>(currentDistance, i));

					pointDistanceMap.erase(pointDistanceMap.rbegin()->first);
				}
			}
			else
			{
				pointDistanceMap.insert(std::pair<double, size_t>(currentDistance, i));
			}
		}
	}

	return pointDistanceMap;
}

// Fill a vector of covariograms over all distances from a given point

std::vector<double> KrigingOperator::CalculateCovariogramVector(const std::vector<SamplePoint>& dataPointCandidate,
	double xCoordinate, double yCoordinate, Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<double> distanceVector(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0);

	double distance = 0.0;
	double covariogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distance = std::sqrt(std::pow(dataPointCandidate[i].x - xCoordinate, 2.0) + std::pow(dataPointCandidate[i].y - yCoordinate, 2.0));

		covariogram = CalculateCovariogram(model, distance, nugget, sill, range);

		distanceVector[i] = covariogram;
	}

	return distanceVector;
}

// Fill a matrix of variograms over all point distances

std::vector<std::vector<double>> KrigingOperator::CalculateVariogramMatrix(const std::vector<SamplePoint>& dataPointCandidate,
	Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<std::vector<double>> distanceMatrix(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), std::vector<double>(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0));

	double distance = 0.0;
	double variogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distanceMatrix[i][i] = CalculateVariogram(model, 0.0, nugget, sill, range);

		for (size_t j = i + 1; j < dataPointCandidate.size(); j++)
		{
			distance = std::sqrt(std::pow(dataPointCandidate[i].x - dataPointCandidate[j].x, 2.0) + std::pow(dataPointCandidate[i].y - dataPointCandidate[j].y, 2.0));

			variogram = CalculateVariogram(model, distance, nugget, sill, range);

			distanceMatrix[i][j] = variogram;
			distanceMatrix[j][i] = distanceMatrix[i][j];
		}
	}

	if (LagrangeMultiplier)
	{
		distanceMatrix[dataPointCandidate.size()][dataPointCandidate.size()] = 0.0;
	}

	return distanceMatrix;
}

// Fill a matrix of covariograms over all point distances

std::vector<std::vector<double>> KrigingOperator::CalculateCovariogramMatrix(const std::vector<SamplePoint>& dataPointCandidate,
	Model model, double nugget, double sill, double range, bool LagrangeMultiplier) const
{
	std::vector<std::vector<double>> distanceMatrix(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), std::vector<double>(LagrangeMultiplier ? dataPointCandidate.size() + 1 : dataPointCandidate.size(), LagrangeMultiplier ? 1.0 : 0.0));

	double distance = 0.0;
	double covariogram = 0.0;

	for (size_t i = 0; i < dataPointCandidate.size(); i++)
	{
		distanceMatrix[i][i] = CalculateCovariogram(model, 0.0, nugget, sill, range);

		for (size_t j = i + 1; j < dataPointCandidate.size(); j++)
		{
			distance = std::sqrt(std::pow(dataPointCandidate[i].x - dataPointCandidate[j].x, 2.0) + std::pow(dataPointCandidate[i].y - dataPointCandidate[j].y, 2.0));

			covariogram = CalculateCovariogram(model, distance, nugget, sill, range);

			distanceMatrix[i][j] = covariogram;
			distanceMatrix[j][i] = distanceMatrix[i][j];
		}
	}

	if (LagrangeMultiplier)
	{
		distanceMatrix[dataPointCandidate.size()][dataPointCandidate.size()] = 0.0;
	}

	return distanceMatrix;
}


void KrigingOperator::CalculatePartialDistanceMap(int start, int end, double* distances, double* variances)
{
	// Rather than take the simple variance of the variable, a better estimation should be the
	// average variance amoung all distances.

	estimatedSill = 0.0;

	double variance = 0.0;
	int count = 0;
	for (size_t i = start; i <= end; i++)
	{
		for (size_t j = i+1; j < _sample_nums; j++)
		{
			variance = std::pow(dataPoint[i].value - dataPoint[j].value, 2.0);

			estimatedSill += variance;
			double distance = std::sqrt(std::pow(dataPoint[i].x - dataPoint[j].x, 2.0) + std::pow(dataPoint[i].y - dataPoint[j].y, 2.0));
			distances[count] = distance;
			variances[count] = variance;
			count++;
		}
	}
    cout<<"rank "<<GetRank()<<" actual pairs: "<<count<<"("<<double(count)/SFi(0,_sample_nums-1)<<")"<<endl;


}

// Fill a map of all distances with semivariogram

void KrigingOperator::CalculateDistanceMap()
{
	// Rather than take the simple variance of the variable, a better estimation should be the
	// average variance amoung all distances.

	estimatedSill = 0.0;

	double variance = 0.0;
	int i0, j0;
	for (size_t i = 0; i < dataPoint.size(); i++)
	{
		for (size_t j = i + 1; j < dataPoint.size(); j++)
		{
			variance = std::pow(dataPoint[i].value - dataPoint[j].value, 2.0);

			estimatedSill += variance;
			double distance = std::sqrt(std::pow(dataPoint[i].x - dataPoint[j].x, 2.0) + std::pow(dataPoint[i].y - dataPoint[j].y, 2.0));
			semiVariogram.insert(std::pair<double, double>(distance, variance));
		}
	}

	estimatedSill = (estimatedSill / 2.0) / ((dataPoint.size() * (dataPoint.size() - 1)) / 2.0);

}


// Find the appropriate lag parameters

void KrigingOperator::CalculateDefaultLagParameters()
{
	// This algorithm seems to come up with reasonable lag parameters.

	double minDistance = semiVariogram.begin()->first;

	lag = minDistance * 2 * 1.25;

	double lagBoundary = (semiVariogram.rbegin()->first - semiVariogram.begin()->first) / 2.0;

	unsigned int lagCount = static_cast<unsigned int>(lagBoundary / lag + 0.5);

	const unsigned MINIMUM_LAG_COUNT = 20;

	if (lagCount < MINIMUM_LAG_COUNT)
	{
		lag = minDistance;
	}

	lagTolerance = static_cast<unsigned int>((lag / 2.0) + 0.5);
}

// Find the estimated variogram parameters

void KrigingOperator::CalculateEstimatedVariogramParameters()
{
	// Note: Determination of sill and range are only rough estimates.
	//       It is up to the user to interpret the empirical and model
	//       variogram plots to assess the validity of the parameters 
	//       used in the chosen model.

	// Rather than take the simple variance of the variable, a better estimation should be the
	// average variance amoung all distances.

	// Sill is the simple variance

	//double mean = std::accumulate(z.begin(), z.end(), 0.0) / z.size();
	//double variance = 0.0;

	//std::for_each(std::begin(z), std::end(z), [&](const double d)
	//{
	//	variance += (d - mean) * (d - mean);
	//});

	//sill = variance / (z.size() - 1);

	// For fixed models, range is the first distance where the sill is reached.

	// For asymptotic models it would be the first distance where the semivariance
	// reaches 95% of the sill.

	// We will assume a fixed model to start.

	estimatedRange = 0.0;

	std::vector<double>::const_iterator semivariance;
	std::vector<double>::const_iterator distance;

	for (semivariance = lagSemivariance.begin(), distance = lagDistance.begin(); semivariance != lagSemivariance.end(); semivariance++, distance++)
	{
		if (*semivariance > estimatedSill)
		{
			estimatedRange = (estimatedRange + *distance) / 2.0;
			break;
		}
		else
		{
			estimatedRange = *distance;
		}
	}

	// Calculate estimated nugget

	double slope = 0.0;	// not used

	DoSimpleLinearRegression(lagDistance, lagSemivariance, &slope, &estimatedNugget);
}

// Simple linear regression

void KrigingOperator::DoSimpleLinearRegression(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double* intercept) const
{
	double Xmean = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
	double Ymean = std::accumulate(Y.begin(), Y.end(), 0.0) / Y.size();

	std::vector<double>::const_iterator Xi = X.begin();
	std::vector<double>::const_iterator Yi = Y.begin();

	double numerator = 0.0;
	double denominator = 0.0;

	while (Xi != lagDistance.end())
	{
		numerator += ((*Xi - Xmean) * (*Yi - Ymean));
		denominator += ((*Xi - Xmean) * (*Xi - Xmean));

		++Xi;
		++Yi;
	}

	*slope = numerator / denominator;
	*intercept = Ymean - (*slope * Xmean);
}

void KrigingOperator::DoSimpleLinearRegressionWithoutIntercept(const std::vector<double>& X, const std::vector<double>& Y, double* slope, double intercept) const
{
	double Xmean = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
	double Ymean = std::accumulate(Y.begin(), Y.end(), 0.0) / Y.size();

	std::vector<double>::const_iterator Xi = X.begin();
	std::vector<double>::const_iterator Yi = Y.begin();

	double numerator = 0.0;
	double denominator = 0.0;

	while (Xi != lagDistance.end())
	{
		numerator += (*Xi * (*Yi - intercept));
		denominator += (*Xi * *Xi);

		++Xi;
		++Yi;
	}

	*slope = numerator / denominator;
}

// Calculate the semivariances

void KrigingOperator::CalculateExperimentalVariogram(double lag, double lagTolerance)
{
	// Clear containers from any previous calculations

	lagDistance.clear();
	lagSemivariance.clear();
	lagCount.clear();

	// Only consider points over half the distance.

	double lagBoundary = (semiVariogram.rbegin()->first - semiVariogram.begin()->first) / 2.0;

	double currentLagDistance = lag / 2.0;

	double currentlagSemivariogram = 0.0;

	unsigned int currentLagCount = 0;

	for (std::multimap<double, double>::const_iterator iterator = semiVariogram.begin(); iterator != semiVariogram.end() && currentLagDistance <= lagBoundary; iterator++)
	{
		if (iterator->first >  currentLagDistance + lagTolerance)
		{
			if (currentLagCount > 0)
			{
				lagDistance.push_back(currentLagDistance);
				lagSemivariance.push_back((currentlagSemivariogram / currentLagCount) / 2.0);
				lagCount.push_back(currentLagCount);

				currentlagSemivariogram = 0.0;
				currentLagCount = 0;
			}

			currentLagDistance += lag;
		}

		if (iterator->first <= currentLagDistance + lagTolerance)
		{
			currentlagSemivariogram += iterator->second;
			currentLagCount++;
		}
	}
}

// Calculate the variogram

double KrigingOperator::CalculateVariogram(Model model, double distance, double nugget, double sill, double range) const
{
	// Notes
	//
	// Linear models do not use nugget, sill or range terminology. For convenience the intercept is passed
	// as the nugget and the intercept as the sill.

	double variogram = 0.0;

	switch (model)
	{
	case Linear:
		variogram = nugget + (sill * distance);
		break;
	case LinearWithoutIntercept:
		variogram = sill * distance;
		break;
	case Spherical:
		if (distance == 0.0)
		{
			variogram = 0.0;
		}
		else if (distance <= range)
		{

			variogram = nugget + (sill - nugget) * ((1.5 * (distance / range) - 0.5 * (std::pow(distance / range, 3.0))));
		}
		else
		{
			variogram = sill;
		}
		break;
	case Exponential:
		if (distance > 0.0)
		{
			variogram = nugget + (sill - nugget) * (1 - std::exp(-(distance / range)));
		}
		else
		{
			variogram = 0.0;
		}
		break;
	case Gaussian:
		if (distance > 0.0)
		{
			variogram = 0.0;
		}
		else
		{
			variogram = nugget + (sill - nugget) * (1 - std::exp(-std::pow(distance / range, 2.0)));
		}
		break;
	default:
		assert(false);
		break;
	}

	return variogram;
}

// Calculate the covariogram

double KrigingOperator::CalculateCovariogram(Model model, double distance, double nugget, double sill, double range) const
{
	// Notes
	//
	// As linear models do not have a sill, it is not possible to calculate a covariogram.

	double covariogram = 0.0;

	switch (model)
	{
	case Linear:
	case LinearWithoutIntercept:
		assert(false);
		break;
	case Spherical:
		if (distance == 0.0)
		{
			covariogram = sill;
		}
		else if (distance <= range)
		{
			covariogram = sill * (1 - (1.5 * (distance / range) - 0.5 * (std::pow(distance / range, 3.0))));
		}
		else
		{
			covariogram = 0.0;
		}
		break;
	case Exponential:
		if (distance == 0.0)
		{
			covariogram = sill;
		}
		else
		{
			covariogram = (sill - nugget) * (std::exp(-distance / range));
		}
		break;
	case Gaussian:
		if (distance == 0.0)
		{
			covariogram = sill;
		}
		else
		{
			covariogram = (sill - nugget) * (std::exp(-std::pow(distance / range, 2.0)));
		}
		break;
	default:
		assert(false);
		break;
	}

	return covariogram;
}

// Calculate models

std::vector<double> KrigingOperator::CalculateLinearModel() const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	// Get slope and intercept

	double slope = 0.0;
	double intercept = 0.0;

	DoSimpleLinearRegression(lagDistance, lagSemivariance, &slope, &intercept);

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Linear, *distance, intercept, slope, 0.0));
	}

	return variogram;
}

std::vector<double> KrigingOperator::CalculateLinearModelWithoutIntercept(double nugget) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	// Get slope

	double slope = 0.0;

	DoSimpleLinearRegressionWithoutIntercept(lagDistance, lagSemivariance, &slope, nugget);

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(LinearWithoutIntercept, *distance, 0.0, slope, 0.0));
	}

	return variogram;
}

std::vector<double> KrigingOperator::CalculateSphericalModel(double nugget, double sill, double range) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Spherical, *distance, nugget, sill, range));
	}

	return variogram;
}

std::vector<double> KrigingOperator::CalculateExponentialModel(double nugget, double sill, double range) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Exponential, *distance, nugget, sill, range));
	}

	return variogram;
}

std::vector<double> KrigingOperator::CalculateGaussianModel(double nugget, double sill, double range) const
{
	std::vector<double> variogram;

	variogram.reserve(lagDistance.size());

	for (std::vector<double>::const_iterator distance = lagDistance.begin(); distance < lagDistance.end(); distance++)
	{
		variogram.push_back(CalculateVariogram(Gaussian, *distance, nugget, sill, range));
	}

	return variogram;
}

// Load the input data


inline double KrigingOperator::getMinDistanceToBlockBound(double x, double y) {
	double boundXY[] = {
		getBlockColIndexByCoord(x)*_blockSize + _sub_extent.minX, //left, block min x
		(getBlockColIndexByCoord(x) + 1)*_blockSize + _sub_extent.minX, //right, block max x
		_sub_extent.maxY - getBlockRowIndexByCoord(y)*_blockSize, //down, block min y
		_sub_extent.maxY - (getBlockRowIndexByCoord(y) + 1)*_blockSize //up, block max y
	};


	double distances[] = {
		abs(x - boundXY[0]),
		abs(x - boundXY[1]),
		abs(y - boundXY[2]),
		abs(y - boundXY[3])
	};
	return min(min(distances[0], distances[1]), min(distances[2], distances[3]));
}

KrigingOperator::~KrigingOperator() {
	vector<SampleBlock>().swap(_pSampleBlocks);
}

int KrigingOperator::readSampleNums(const char* filename, char** pSpatialRefWkt) {
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
	OGRSpatialReference* sref = poLayer->GetSpatialRef();
	sref->exportToWkt(pSpatialRefWkt);

	OGRFeature* poFeature;

	poLayer->ResetReading();
	while ((poFeature = poLayer->GetNextFeature()) != NULL) {
		_sample_nums++;
		OGRFeature::DestroyFeature(poFeature);
	}
	//_sample_nums = poLayer->GetFeatureCount();//why not use this
#if GDAL_VERSION_MAJOR >= 2
	GDALClose(poDatasetsrc);
#else
	OGRDataSource::DestroyDataSource(poDatasetsrc);
#endif
	return _sample_nums;
}

bool KrigingOperator::readSamples(const char* filename, int fieldIdx, char** pSpatialRefWkt, vector<SamplePoint> &samples) {

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
		OGRGeometry* poGeometry;
		poGeometry = poFeature->GetGeometryRef();
		if (poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint) {
			//why not poLayer->GetExtent
			OGRPoint* poPoint = (OGRPoint *)poGeometry;
			x = poPoint->getX();
			y = poPoint->getY();
			SamplePoint point;
			point.x = x;
			point.y = y;
			point.value = poFeature->GetFieldAsDouble(fieldIdx);
			//OGRFieldDefn* name = poFeature->GetFieldDefnRef(fieldIdx);
			//cout << name->GetNameRef() << endl;
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


void KrigingOperator::creatSampleBlocks(vector<SamplePoint> &samples) {
	dataPoint = samples;
	_blockRows = ceil((_glb_extent.maxY - _glb_extent.minY) / _blockSize);
	_blockCols = ceil((_glb_extent.maxX - _glb_extent.minX) / _blockSize);
	_pSampleBlocks.resize(_blockRows * _blockCols);
	for (int i = 0; i < _sample_nums; ++i) {
		double x = samples[i].x;
		double y = samples[i].y;
		int iCol = getBlockColIndexByCoord(x);
		int iRow = getBlockRowIndexByCoord(y);
		int index = iRow * _blockCols + iCol;
		if (index >= _pSampleBlocks.size()) {
			cerr << "Index out of range. An error in creatSampleBlocks()" << endl;
			cout << "__blockRows * _blockCols = " << _blockRows << " * " << _blockCols << endl;
			printf("_pSampleBlocks.size()=%llu,samples.size()=%llu, _sample_nums=%i\n", _pSampleBlocks.size(), samples.size(), _sample_nums);
		}
		_pSampleBlocks[index].samplePoints.push_back(samples[i]);
	}
}
void KrigingOperator::idwLayer(RasterLayer<double>& layerD, char** pSpatialRefWkt, DomDcmpType dcmpType) {
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
	layerD._pMetaData->_domDcmpType = ROWWISE_DCMP;
	MPI_Comm_rank(MPI_COMM_WORLD, &layerD._pMetaData->myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &layerD._pMetaData->processor_number);
	_myRank = layerD._pMetaData->myrank;

	DeComposition<double> deComp(layerD._pMetaData->_glbDims, *(layerD.nbrhood()));

	deComp.rowDcmp(*(layerD._pMetaData), layerD._pMetaData->processor_number);

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
void KrigingOperator::maskLayer(RasterLayer<int>& layerD) {
	_pMaskLayer = &layerD;
}
void KrigingOperator::initIdwLayerGlobalInfo(RasterLayer<double>& layerD, char** pSpatialRefWkt) {
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
	layerD._pMetaData->_domDcmpType = ROWWISE_DCMP;
	MPI_Comm_rank(MPI_COMM_WORLD, &layerD._pMetaData->myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &layerD._pMetaData->processor_number);
	_myRank = layerD._pMetaData->myrank;

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

}
void KrigingOperator::idwLayer(RasterLayer<double>& layerD, char** pSpatialRefWkt, CoordBR &subWorkBR) {

	initIdwLayerGlobalInfo(layerD, pSpatialRefWkt);


	layerD._pMetaData->_MBR = subWorkBR;
	layerD._pMetaData->_localdims = SpaceDims(subWorkBR.nRows(), subWorkBR.nCols());
	layerD.nbrhood()->calcWorkBR(layerD._pMetaData->_localworkBR, layerD._pMetaData->_localdims);
	layerD.newCellSpace(layerD._pMetaData->_localdims); //read global samples, only calculate local workBR

	_sub_extent.minX = _glb_extent.minX;
	_sub_extent.maxX = _glb_extent.maxX;
	_sub_extent.maxY = _glb_extent.maxY - layerD._pMetaData->_MBR.minIRow() * _cellSize;
	_sub_extent.minY = _glb_extent.maxY - layerD._pMetaData->_MBR.maxIRow() * _cellSize - _cellSize;
	_xSize = layerD._pMetaData->_localdims.nCols();
	_ySize = layerD._pMetaData->_localdims.nRows();

	_pIDWLayer = &layerD;
	Configure(_pIDWLayer, false);
}
void KrigingOperator::idwLayerSerial(RasterLayer<double>& layerD, char** pSpatialRefWkt) {
	if (GetRank() != 0) {
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
bool KrigingOperator::isTermination() {
	return flag;
}

std::multimap<double, SamplePoint> KrigingOperator::searchNbrSamples(const int subMinRow, int cellRow, int cellCol, double cellX, double cellY) {
	int blockRow = getBlockRowIndexByCoord(cellY);
	int blockRow1 = getBlockRowIndexByCellIndex(cellRow + subMinRow);
	int blockCol = getBlockColIndexByCoord(cellX);
	int blockRows = _blockRows;
	int blockCols = _blockCols;

	std::multimap<double, SamplePoint> candidateDistanceMap;
	int searchRad = 0; //search circle by circle. 1 = 3*3 neighbor
	bool isSearch = true;
	double minDistToBound = getMinDistanceToBlockBound(cellX, cellY);
	while (isSearch) {
		double searchRange = (double)(searchRad - 1) * _blockSize + minDistToBound; // actual search range
		if (_idw_buffer>0 && searchRange >= _idw_buffer) {
			break;
		}
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
				//store all of the beginning and ending rows
				for (int tCol = blockCol - searchRad; tCol <= blockCol + searchRad; ++tCol) {
					if (tCol < 0 || tCol >= blockCols) {
						continue;
					}
					block2search[blockCount++] = tRow * blockCols + tCol;
				}
			}
			else {
				if (blockCol - searchRad >= 0 && blockCol - searchRad < blockCols) {
					block2search[blockCount++] = tRow * blockCols + blockCol - searchRad;
				}
				if (blockCol + searchRad >= 0 && blockCol + searchRad < blockCols) {
					block2search[blockCount++] = tRow * blockCols + blockCol + searchRad;
				}
			}
		}

		//traverse int block2search[2*searchRad*4] block idx;
		//search _pSampleBlocks[i].samplePoints
		for (int i = 0; i < blockCount; ++i) {

			int blockIdx = block2search[i];

			for (vector<SamplePoint>::iterator iter = _pSampleBlocks[blockIdx].samplePoints.begin(); iter != _pSampleBlocks[blockIdx].samplePoints.end(); ++iter) {

				double currentDistance = sqrt(pow(iter->x - cellX, 2) + pow(iter->y - cellY, 2));

				if (candidateDistanceMap.size() == _nbrPoints)
				{
					if (currentDistance < candidateDistanceMap.rbegin()->first)
					{
						candidateDistanceMap.insert(std::pair<double,SamplePoint>(currentDistance, *iter));

						candidateDistanceMap.erase(candidateDistanceMap.rbegin()->first);
					}
				}
				else
				{
					candidateDistanceMap.insert(std::pair<double, SamplePoint>(currentDistance, *iter));
				}


			}
		}

		if (candidateDistanceMap.size() == _nbrPoints) {
            //points are enough, and actual search range > max distance among point candidates
			if (searchRange >= candidateDistanceMap.rbegin()->first) {
				isSearch = false;
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
    
	return candidateDistanceMap;
}


bool KrigingOperator::Operator(const CellCoord& coord, bool operFlag) {
	double startTime = MPI_Wtime();
	int iRow = coord.iRow();
	int iCol = coord.iCol();

	double mask = 1;
	double maskNoData = 0;
	if (_pMaskLayer) {
		int maskRow = _pIDWLayer->rowAtOtherLayer(_pMaskLayer, iRow);
		int maskCol = _pIDWLayer->colAtOtherLayer(_pMaskLayer, iCol);
        if(maskRow>=_pMaskLayer->cellSpace()->nRows()) {
            if(iCol==0) {
                cout<<"Reading mask layer row "<<maskRow<<" exceeding boundary "<<_pMaskLayer->cellSpace()->nRows()<<", converted to boundary"<<endl;
                double scale = (double)_pMaskLayer->metaData()->_glbDims.nRows()/_pIDWLayer->metaData()->_glbDims.nRows();
                cout<< scale*iRow<<endl;
            }
            maskRow=_pMaskLayer->cellSpace()->nRows()-1;
        }
        if(maskCol>=_pMaskLayer->cellSpace()->nCols()) {
            if(iRow==0) {
                cout<<"Reading mask layer col "<<maskCol<<" exceeding boundary "<<_pMaskLayer->cellSpace()->nCols()<<", converted to boundary"<<endl;
            }
            maskCol=_pMaskLayer->cellSpace()->nCols()-1;
        }
		mask = (*_pMaskLayer->cellSpace())[maskRow][maskCol];
		maskNoData = _pMaskLayer->metaData()->noData;
	}

	if (mask == maskNoData) {
		if (_pComptLayer) {
			(*_pIDWLayer->cellSpace())[iRow][iCol] = (MPI_Wtime() - startTime) * 1000;
		}else {
		    (*_pIDWLayer->cellSpace())[iRow][iCol] = _noData;
		}
		return true;
	}

	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	CellSpace<double>& idwL = *_pIDWLayer->cellSpace();
	const int minRow = _pIDWLayer->_pMetaData->_MBR.minIRow();


	double cellX = getXByCellIndex(iCol);
	double cellY = getYByCellIndex(iRow + minRow);
	multimap<double, SamplePoint> candidateDistanceMap = searchNbrSamples(minRow, iRow, iCol, cellX, cellY);

	vector<SamplePoint> dataPointCandidate;
	for (std::multimap<double, SamplePoint>::const_iterator iterator = candidateDistanceMap.begin();iterator != candidateDistanceMap.end();iterator++)
	{
		dataPointCandidate.push_back(iterator->second);
	}
	vector<double> distanceCovariogramVector = CalculateCovariogramVector(dataPointCandidate, cellX, cellY, KrigingOperator::Spherical, GetEstimatedNugget(), GetEstimatedSill(), GetEstimatedRange(), true);

	vector<vector<double>> distanceCovariogramMatrix = CalculateCovariogramMatrix(dataPointCandidate, KrigingOperator::Spherical, GetEstimatedNugget(), GetEstimatedSill(), GetEstimatedRange(), true);

	LUDecomposition luDecomposition(distanceCovariogramMatrix);
	luDecomposition.Decompose();

	vector<double> weights = luDecomposition.Solve(distanceCovariogramVector);
	double estimatedZ = std::inner_product(weights.begin(), weights.end() - 1,
		dataPointCandidate.begin(),
		0.0, plus<double>(),
		[](double weight, const SamplePoint& dataPointCandidate) { return weight * dataPointCandidate.value; });
	//cout << estimatedZ; // -nan(ind)
	idwL[iRow][iCol] = estimatedZ;
	double endTime = MPI_Wtime();
	if (_pComptLayer) {
		idwL[iRow][iCol] = (endTime - startTime) * 1000;
	}
	return true;
}
