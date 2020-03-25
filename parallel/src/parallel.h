#pragma once

#define ROOTPROC 0
#include <math.h>
#include <mpi.h>
#include "matrix_helpers.h"
#include "config.h"


void premultA(int row_size, int row_beginning, int K, int Ny, double dx2, double dy2, double *r0_l, double *x);

// Very specific parallel gradient descent just for this problem

void pcgd(  int rank, int size,         ///< MPI formulation
            double dx2, double dy2,     ///< Problem formulation
            int Ny, int K,              ///< A/b matrix formulation
            double* x, double* b);      ///< A/b problem     