#ifndef Application_H
#define Application_H

/***************************************************************************
* application.h
*
* Project: GPRO, v 1.0
* Purpose: Header file for class GPRO::Application
* Author:  Zhan Lijun
* E-mail:  zhanlj@lreis.ac.cn
****************************************************************************
* Copyright (c) 2013. Zhan Lijun
* NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
* purposes, NO COMMERCIAL usages are allowed unless the author is 
* contacted and a permission is granted
* 
****************************************************************************/

#include "basicTypes.h"
#include "mpi.h"

namespace GPRO
{
  class Application 
  {
    public:
	  static bool START(ProgramType programType, int argc, char *argv[]);
      static bool END();

	public:
      static ProgramType _programType;
  };
};


#endif