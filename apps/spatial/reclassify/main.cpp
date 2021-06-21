#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computeLayer.h"
#include "application.h"
#include "reclassifyOperator.h"
#include "transformation.h"
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
    char* dataNeighbor;
    char* outputFileName;
    int nodataLoad;
    int validLoad;
    bool decomposeBySapce; /// decomp by compute load if false
    char* writeLoadPath=nullptr;

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
                dataNeighbor = argv[i];
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
        if (strcmp(argv[i], "-nodataLoad") == 0) {
            i++;
            if (argc > i) {
                nodataLoad = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-nodataLoad'!");
            }
        }
        if (strcmp(argv[i], "-validLoad") == 0) {
            i++;
            if (argc > i) {
                validLoad = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-validLoad'!");
            }
        }
        if (strcmp(argv[i], "-dcmp") == 0) {
            i++;
            if (argc > i) {
                decomposeBySapce = ifDecomposeBySpace(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-dcmp'!");
            }
        }
        if (strcmp(argv[i], "-writeLoad") == 0) {
            i++;
            if (argc > i) {
                writeLoadPath = argv[i];
                i++;
            }
            else {
                Usage("No argument followed '-writeLoad'!");
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

        int name_len = MPI_MAX_PROCESSOR_NAME;
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        MPI_Get_processor_name(processor_name,&name_len);
        if (myRank == 0) {
            cout<<"PaRGO-Reclassify. "<<process_nums<<" core(s)"<<endl;
            for (int i = 0; i < argc; ++i) {
                cout<<argv[i];
                if(argv[i][0]=='-') 
                    cout<<" ";
                else
                    cout<<endl;
            }
            cout<<endl;
        }

        double starttime = MPI_Wtime();
        double endtime;
        
        if (decomposeBySapce) {
            RasterLayer<double> inputLayer("inputLayer");
            inputLayer.readFile(inputFileName, ROWWISE_DCMP);

            RasterLayer<double> outputLayer("outputLayer");
            outputLayer.copyLayerInfo(inputLayer);

            ReclassifyOperator recOper;
            recOper.setInputLayer(inputLayer);
            recOper.setOutputLayer(outputLayer);
            recOper.Run();
            outputLayer.writeFile(outputFileName);
        } else {
            // ONLY support write load by intensity ratio
             vector<RasterLayer<double>*> inputLayers;
             RasterLayer<double> inputLayerGlobal;
             inputLayerGlobal.readGlobalFileSerial(inputFileName);
             inputLayers.emplace_back(&inputLayerGlobal);
             CoordBR subWorkBR;
             starttime = MPI_Wtime();
             ComputeLayer<double> comptLayer("computLayer");
             comptLayer.initSerial(inputLayers,dataNeighbor,1);
             Transformation<double> transOper(nodataLoad, validLoad, &comptLayer); 
             transOper.run();
             if (writeLoadPath) {
                 comptLayer.writeComputeIntensityFileSerial(writeLoadPath);               
             }
             comptLayer.getCompuLoad(ROWWISE_DCMP, process_nums, subWorkBR);
             if (myRank == 0)
                 cout << "compute-write-dcmp time is " << MPI_Wtime() - starttime << endl;


        }
	    cout<<"run time is "<<MPI_Wtime()-starttime<<endl;
        Application::END();
        return 0;
    }
}