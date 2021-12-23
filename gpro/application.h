/**
 * \file application
 * \author Zhan Lijun (zhanlj@lreis.ac.cn)
 * \brief Header file for class GPRO::Application
 * \version 1.0
 * 
 * \copyright Copyright (c) 2013-2020
 *  NOTE: this library can ONLY be used for EDUCATIONAL and SCIENTIFIC 
 *  purposes, NO COMMERCIAL usages are allowed unless the author is 
 *  contacted and a permission is granted
 * 
 * changelog:
 *  - 1. 2020 - Wang Yujing - Code reformat
 */

#ifndef Application_H
#define Application_H

#include "basicTypes.h"
#include "mpi.h"

namespace GPRO {
    /**
     * \ingroup gpro
     * \class Application
     * \brief To initialize or finalize outside frameworks such as MPI 
     */
    class Application {
    public:
        /**
         * \brief initialize outside framework
         * \param[in] programType type of outside framework
         * \param[in] argc same as the main function
         * \param[in] argv same as the main function
         */
        static bool START(ProgramType programType, int argc, char* argv[]);
        /**
        * \brief finalize outside framework
        */
        static bool END();

    public:
        static ProgramType _programType; //serial, MPI, openMP or CUDA program
    };
};

#endif
