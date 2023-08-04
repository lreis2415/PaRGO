#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computeLayer.h"
#include "application.h"
#include "krigingOperator.h"
#include "krigingTransformation.h"
using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: kriging -input <input raster file> -output <reclassified file> " << endl;

    exit(1);
}

int main(int argc, char* argv[]) {

	char* inputFileName=nullptr;
	char* maskFileName = nullptr;
	char* outputFileName = nullptr;
	char* dataNeighbor = nullptr;
	char* compuNeighbor = nullptr;
	int fieldIndex;
	float cellSize;
	int searchPoints;
	double searchRange;
	double blockSize; //unit: meter
	int granularity = 10;//resolution of the computational domain = granularity * resolution of the data domain (10 by default)
	bool decomposeBySapce; /// decomp by compute load if false
	char* writeLoadPath = nullptr;
	char* readLoadPath = nullptr;

	int i = 1;
	bool simpleusage = true;
	while (argc > i) {
		if (strcmp(argv[i], "-sample") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				inputFileName = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-sample'!");
			}
		}
		if (strcmp(argv[i], "-mask") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				maskFileName = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-mask'!");
			}
		}
		else if (strcmp(argv[i], "-dataNbr") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				dataNeighbor = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-dataNbr'!");
			}
		}
		else if (strcmp(argv[i], "-computeNbr") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				compuNeighbor = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-computeNbr'!");
			}
		}
		else if (strcmp(argv[i], "-out") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				outputFileName = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-out'!");
			}
		}
		else if (strcmp(argv[i], "-resolution") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				cellSize = atof(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-resolution'!");
			}
		}
		else if (strcmp(argv[i], "-fieldIndex") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				fieldIndex = atoi(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-fieldIndex'!");
			}
		}
		else if (strcmp(argv[i], "-searchPointNum") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				searchPoints = atoi(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-searchPointNum'!");
			}
		}
		else if (strcmp(argv[i], "-searchRange") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				searchRange = atof(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-searchRange'!");
			}
		}
		else if (strcmp(argv[i], "-blockSize") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				blockSize = atof(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-blockSize'!");
			}
		}
		else if (strcmp(argv[i], "-granularity") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				granularity = atoi(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-granularity'!");
			}
		}
		else if (strcmp(argv[i], "-dcmp") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				decomposeBySapce = ifDecomposeBySpace(argv[i]);
				i++;
			}
			else {
				Usage("No argument followed '-dcmp'!");
			}
		}
		else if (strcmp(argv[i], "-writeLoad") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				writeLoadPath = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-writeLoad'!");
			}
		}
		else if (strcmp(argv[i], "-readLoad") == 0) {
			simpleusage = false;
			i++;
			if (argc > i) {
				readLoadPath = argv[i];
				i++;
			}
			else {
				Usage("No argument followed '-readLoad'!");
			}
		}
		else {
			// Simple Usage
			if (!simpleusage) Usage();
		}
	}

	Application::START(MPI_Type, argc, argv);
	int myRank, process_nums;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
	int name_len = MPI_MAX_PROCESSOR_NAME;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(processor_name, &name_len);
	if (myRank == 0) {
		cout << "PaRGO-Kriging. " << process_nums << " core(s)" << endl;
		for (int i = 0; i < argc; ++i) {
			cout << argv[i];
			if (argv[i][0] == '-')
				cout << " ";
			else
				cout << endl;
		}
		cout << endl;
	}
	// cout << "process " << myRank << " on " << processor_name << endl;

	double starttime;


	KrigingOperator krigingOper(cellSize, searchPoints, searchRange, blockSize);
	char* spatialrefWkt;

	krigingOper.readSampleNums(inputFileName, &spatialrefWkt);
	vector<SamplePoint> samples;
	krigingOper.readSamples(inputFileName, fieldIndex, &spatialrefWkt, samples);
	krigingOper.creatSampleBlocks(samples);

	starttime = MPI_Wtime();
	krigingOper.Initialize();
	cout << "kriging serial init time is " << MPI_Wtime() - starttime << endl << endl;


	// starttime = MPI_Wtime();
	// krigingOper.ParallelInitialize();
	// cout << "kriging parallel init time is " << MPI_Wtime() - starttime << endl << endl;


    //不删除数组，存到类里后面要用
	//vector<SamplePoint> vTemp(0);
	//vTemp.swap(samples);
	RasterLayer<int> maskLayer("maskLayer");
	if (maskFileName!=nullptr) {
		maskLayer.readNeighborhood(dataNeighbor);
	}
	RasterLayer<double> idwLayer("idwLayer");
	idwLayer.readNeighborhood(dataNeighbor);

	if (decomposeBySapce) {
		krigingOper.idwLayer(idwLayer, &spatialrefWkt);
		if (maskFileName != nullptr) {
			maskLayer.readFile(maskFileName, ROWWISE_DCMP);
			krigingOper.maskLayer(maskLayer);
		}
		ComputeLayer<double> comptLayer("computeLayer");
		//if (writeLoadPath) {
		//	comptLayer.addRasterLayerSerial(&idwLayer);
		//	comptLayer.initSerial(compuNeighbor, 1);
		//	krigingOper.comptLayer(comptLayer);
		//}

        // Actually not in use! The output of parallel Preliminary-Experiment-Computational-Domain-Constructing is in the "-out" path.

        if (writeLoadPath) {
			comptLayer.addRasterLayer(&idwLayer);
			comptLayer.init(compuNeighbor, 1);
			krigingOper.comptLayer(comptLayer);
		}
		if (myRank == 0) cout << "start computing" << endl;

		starttime = MPI_Wtime();
		krigingOper.Run(); //fill idwLayer
		//if (writeLoadPath)
		//	comptLayer.writeComputeIntensityFileSerial(writeLoadPath);
      //  if (writeLoadPath)
		    //comptLayer.writeComputeIntensityFile(writeLoadPath);
	}
	else {
		krigingOper.idwLayerSerial(idwLayer, &spatialrefWkt);
		CoordBR subWorkBR;

		if (myRank == 0) cout << "start dcmp" << endl;
		if (readLoadPath) {
			starttime = MPI_Wtime();
			ComputeLayer<double> comptLayer("computLayer");
			comptLayer.addRasterLayerSerial(&idwLayer);
			comptLayer.init(compuNeighbor, granularity);
			comptLayer.readComputeLoadFile(readLoadPath);
			comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
		    if (myRank == 0) cout << "dcmp-read time is " << MPI_Wtime() - starttime << endl;
		}
		else {
			RasterLayer<int> fullMaskLayer("fullMaskLayer");
			fullMaskLayer.readNeighborhood(dataNeighbor);
			fullMaskLayer.readGlobalFileSerial(maskFileName);
			krigingOper.maskLayer(fullMaskLayer);
			starttime = MPI_Wtime();
			ComputeLayer<double> comptLayer("computLayer");
			comptLayer.addRasterLayerSerial(&idwLayer);
			comptLayer.init(compuNeighbor, granularity);
			KrigingTransformation trans(&comptLayer, &krigingOper);
			trans.run();
			comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
		    if (myRank == 0) cout << "dcmp-write time is " << MPI_Wtime() - starttime << endl;
			if (writeLoadPath) {
				comptLayer.writeComputeIntensityFile(writeLoadPath);
			}
		}
		cout << myRank << " subWorkBR " << subWorkBR.minIRow() << " " << subWorkBR.maxIRow() << " " << subWorkBR.nRows() << endl;
		krigingOper.idwLayer(idwLayer, &spatialrefWkt, subWorkBR);
		maskLayer.readFile(maskFileName, subWorkBR, ROWWISE_DCMP);
		krigingOper.maskLayer(maskLayer);

		if (myRank == 0) cout << "start computing" << endl;
		starttime = MPI_Wtime();
		krigingOper.Run();
	}
	if (myRank == 0)
		cout << "compute time is " << MPI_Wtime() - starttime << endl << endl;
	idwLayer.writeFile(outputFileName);

	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "write done." << endl;

	Application::END();
	return 0;
}