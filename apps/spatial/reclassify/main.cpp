#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computeLayer.h"
#include "application.h"
#include "reclassifyOperator.h"
using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << " Usage: reclassify -input <input raster file> -output <reclassified file> " << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName;
    char* outputFileName;

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
        if (strcmp(argv[i], "-nbr") == 0) {
            i++;
            if (argc > i) {
                inputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-nbr'!");
            }
        }
        if (strcmp(argv[i], "-output") == 0) {
            i++;
            if (argc > i) {
                outputFileName = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-output'!");
            }
        }
        if (!FileExists(inputFileName)) {
            Usage("input file not exists");
        }
        if (!outputFileName) {
            Usage("please specify an output file");
        }
        Application::START(MPI_Type, argc, argv);

        double starttime=MPI_Wtime();

        RasterLayer<double> inputLayer("inputLayer");
        inputLayer.readFile(inputFileName, ROWWISE_DCMP);

        RasterLayer<double> outputLayer("outputLayer");
        outputLayer.copyLayerInfo(inputLayer);

        ReclassifyOperator recOper;
        recOper.setInputLayer(inputLayer);
        recOper.setOutputLayer(outputLayer);
        recOper.Run();
        outputLayer.writeFile(outputFileName);
	    cout<<"run time is "<<MPI_Wtime()-starttime<<endl;
        Application::END();
        return 0;
    }
}