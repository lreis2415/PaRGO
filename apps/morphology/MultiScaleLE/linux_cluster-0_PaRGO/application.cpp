#include"application.h"

GPRO::ProgramType GPRO::Application::_programType = Serial_Type;


bool GPRO::Application::
START(ProgramType programType, int argc, char *argv[])
{
	_programType = programType;

	if(_programType == MPI_Type || _programType == MPI_OpenMP_Type)
	{
		MPI_Init(&argc, &argv);
	}
	else if(_programType == CUDA_Type)
	{
		;
	}
	else
	{
		;
	}
	
	
	return true;
}

bool GPRO::Application::
END()
{
	if(_programType == MPI_Type || _programType == MPI_OpenMP_Type)
	{
		MPI_Finalize();
	}
	else if(_programType == CUDA_Type)
	{
		;
	}
	else
	{
		;
	}
	
}