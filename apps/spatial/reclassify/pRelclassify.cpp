#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computeLayer.h"
#include "application.h"
#include "reclassifyOperator.h"
#include "transformation.h"
#include <algorithm>
#include <map>
using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: reclassify -input <input raster file> -output <reclassified file> " << endl;
    cout << "Example.1. reclassify -input /path/to/input.tif -output /path/to/output.tif -levels 100,150,200,300" << endl;
    cout << "Example.1. reclassify -input /path/to/input.tif -output /path/to/moore.nbr -levels 100,150,200,300 -dcmp compute -nodataLoad 0 -validLoad 1" << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName = nullptr;
    char* outputFileName = nullptr;
    int nodataLoad;
    int validLoad;
    bool decomposeBySapce = true; /// decomp by compute load if false
    char* writeLoadPath = nullptr;
    char* cLevels = nullptr;

    int i = 1;
    while (argc > i) {
        if (strcmp(argv[i], "-input") == 0) {
            i++;
            if (argc > i) {
                inputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-input'!");
            }
        }
        else if (strcmp(argv[i], "-output") == 0) {
            i++;
            if (argc > i) {
                outputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-output'!");
            }
        }
        else if (strcmp(argv[i], "-levels") == 0) {
            i++;
            if (argc > i) {
                cLevels = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-levels'!");
            }
        }
        else if (strcmp(argv[i], "-nodataLoad") == 0) {
            i++;
            if (argc > i) {
                nodataLoad = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-nodataLoad'!");
            }
        }
        else if (strcmp(argv[i], "-validLoad") == 0) {
            i++;
            if (argc > i) {
                validLoad = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-validLoad'!");
            }
        }
        else if (strcmp(argv[i], "-dcmp") == 0) {
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
            i++;
            if (argc > i) {
                writeLoadPath = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-writeLoad'!");
            }
        }
    }
    if (!FileExists(inputFileName)) {
        Usage("input file not exists");
    }
    if (!outputFileName) {
        Usage("please specify an output file");
    }
    Application::START(MPI_Type, argc, argv);

    int myRank, process_nums;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_nums);
    
#ifdef _DEBUG
    int name_len = MPI_MAX_PROCESSOR_NAME;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);
    if (myRank == 0) {
        cout << "PaRGO-Reclassify. " << process_nums << " core(s) by "<<processor_name << endl;
        for (int i = 0; i < argc; ++i) {
            cout << argv[i];
            if (argv[i][0] == '-')
                cout << " ";
            else
                cout << endl;
        }
        cout << endl;
    }
#endif

    vector<string> tokens = SplitString(cLevels, ',');
    vector<double> dLevels;
    for (int i = 0; i < tokens.size(); ++i) {
        dLevels.emplace_back(stod(tokens[i]));
    }
    sort(dLevels.begin(), dLevels.end());
    

    double starttime = MPI_Wtime();
    double endtime;
    // no load-balancing, equal-area strategy
    if (decomposeBySapce) {
        RasterLayer<double> inputLayer("inputLayer");
        inputLayer.readFile(inputFileName, ROWWISE_DCMP);

        RasterLayer<double> outputLayer("outputLayer");
        outputLayer.copyLayerInfo(inputLayer);

        ReclassifyOperator recOper;
        recOper.setInputLayer(inputLayer);
        recOper.setOutputLayer(outputLayer);
        recOper.setLevels(&dLevels);
        recOper.Run();
        outputLayer.writeFile(outputFileName);
    }
    else {
        // load-balancing, intensity ratio mode
        // ONLY support write load by intensity ratio
        RasterLayer<double> inputLayerGlobal;
        inputLayerGlobal.readGlobalFileSerial(inputFileName);
        starttime = MPI_Wtime();
        ComputeLayer<double> comptLayer("computLayer");
        comptLayer.init(&inputLayerGlobal, nullptr, 10);
        Transformation<double> transOper(nodataLoad, validLoad, &comptLayer);
        transOper.run();
        if (writeLoadPath) {
            comptLayer.writeComputeIntensityFile(writeLoadPath);
        }
        CoordBR subWorkBR;
        comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
        if (myRank == 0)
            cout << "compute-write-dcmp time is " << MPI_Wtime() - starttime << endl;

        RasterLayer<double> inputLayer("inputLayer");
        inputLayer.readFile(inputFileName, subWorkBR);

        RasterLayer<double> outputLayer("outputLayer");
        outputLayer.copyLayerInfo(inputLayer);

        ReclassifyOperator recOper;
        recOper.setInputLayer(inputLayer);
        recOper.setOutputLayer(outputLayer);
        recOper.setLevels(&dLevels);
        recOper.Run();
        outputLayer.writeFile(outputFileName);

    }
    cout << "compute time is " << MPI_Wtime() - starttime << endl;
    Application::END();
    return 0;
}
