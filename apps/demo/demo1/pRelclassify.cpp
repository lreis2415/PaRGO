#include <iostream>
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "reclassifyOperator.h"
#include "transformation.h"
using namespace std;
using namespace GPRO;

/*!
 * The main function is the entrance of your algorithm.
 * You can read the files, parameters, and run your Operator here.
 */

int main(int argc, char* argv[]) {
    /*!
	 * Parse input arguments.
	 * DO NOT start the application unless the required inputs are provided!
	 */

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
        }
        else if (strcmp(argv[i], "-output") == 0) {
            i++;
            if (argc > i) {
                outputFileName = argv[i];
                i++;
            }
        }
    }

    Application::START(MPI_Type, argc, argv); //Start the PaRGO application.

    RasterLayer<double> inputLayer("inputLayer"); //declare the input layer
    inputLayer.readFile(inputFileName, ROWWISE_DCMP); //read input data using row-wise decomposition

    RasterLayer<double> outputLayer("outputLayer"); //declare the output layer
    outputLayer.copyLayerInfo(inputLayer); //initialize the output layer using the metadata of input layer

    vector<double> dLevels; //set the reclassify values
    dLevels.emplace_back(-100);
    dLevels.emplace_back(50);
    dLevels.emplace_back(100);
    dLevels.emplace_back(200);
    dLevels.emplace_back(300);
    dLevels.emplace_back(9999);

    ReclassifyOperator recOper; //declare the operator.
    recOper.setInputLayer(inputLayer); //pass the input RasterLayer to the operator.
    recOper.setOutputLayer(outputLayer); //pass the output RasterLayer to the operator.
    recOper.setLevels(&dLevels); //pass other parameters the operator needs.

    recOper.Run(); //run the operator.
    
    outputLayer.writeFile(outputFileName); // write output.
    cout<<"write done."<<endl;

    Application::END(); //finish the PaRGO application. Release resources.

    return 0;
}
