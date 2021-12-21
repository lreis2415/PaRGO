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


int main(int argc, char* argv[]) {
    char* inputFileName = nullptr;
    char* outputFileName = nullptr;

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
    }
    if (!FileExists(inputFileName)) {
        Usage("input file does not exist");
    }
    if (!outputFileName) {
        Usage("please specify an output file");
    }

    Application::START(MPI_Type, argc, argv);

    RasterLayer<double> inputLayer("inputLayer");
    inputLayer.readFile(inputFileName, ROWWISE_DCMP);

    RasterLayer<double> outputLayer("outputLayer");
    outputLayer.copyLayerInfo(inputLayer);

    vector<double> dLevels={50,100,200,300};

    ReclassifyOperator recOper;
    recOper.setInputLayer(inputLayer);
    recOper.setOutputLayer(outputLayer);
    recOper.setLevels(&dLevels);

    recOper.Run();
    
    outputLayer.writeFile(outputFileName);

    Application::END();
    return 0;
}
