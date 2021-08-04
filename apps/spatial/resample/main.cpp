#include <iostream>
#include <string>
#include "mpi.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "computeLayer.h"
#include "application.h"
#include "resampleOperator.h"
#include "transformation.h"
using namespace std;
using namespace GPRO;

void Usage(const string& error_msg = "") {
    if (!error_msg.empty()) {
        cout << "FAILURE: " << error_msg << endl << endl;
    }

    cout << "NOTE: This operator is for preliminary experiment mode. Only support upscaling by a specified granularity."<<endl
    <<"e.g., -g 2 means one output cell is determined by 2*2 cells in the input (resolution from 10m to 20m). "<<endl<<endl
    <<"Usage: resample -input <input raster file> -output <resampled file> -g <granularity, 10 if not specified>" << endl
    <<"TODO: the output resolution (cellsize) in metadata is not changed yet." << endl;

    exit(1);
}

int main(int argc, char* argv[]) {
    char* inputFileName;
    char* outputFileName;
    int g;

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
        if (strcmp(argv[i], "-g") == 0) {
            i++;
            if (argc > i) {
                g = atoi(argv[i]);
                i++;
            }
            else {
                Usage("No argument followed '-g'!");
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
            cout<<"PaRGO-Resample. "<<process_nums<<" core(s)"<<endl;
            for (int i = 0; i < argc; ++i) {
                cout<<argv[i];
                if(argv[i][0]=='-') 
                    cout<<" ";
                else
                    cout<<endl;
            }
            cout<<endl;
        }

        double endtime;
        
        RasterLayer<double> inputLayer("inputLayer");
        inputLayer.newRectangleNbrhood(3*g);
        inputLayer.readFile(inputFileName, ROWWISE_DCMP);

        RasterLayer<double> outputLayer("outputLayer");
        outputLayer.newLocalNbrhood();
        outputLayer.newUpscaleFile(inputFileName,g,ROWWISE_DCMP);

        ResampleOperator reOper;
        reOper.setInputLayer(inputLayer);
        reOper.setOutputLayer(outputLayer);
        reOper.setG(g);
        double starttime = MPI_Wtime();
        reOper.Run();
        MPI_Barrier(MPI_COMM_WORLD);
        if(myRank==0) {
	        cout<<"run time is "<<MPI_Wtime()-starttime<<endl;
        }
        outputLayer.writeFile(outputFileName);

        

        Application::END();
        return 0;
    }
}