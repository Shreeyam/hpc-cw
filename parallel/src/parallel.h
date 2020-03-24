#pragma once

#include <thread>
#include <chrono>
#include <math.h>
#include <mpi.h>
#include "matrix_helpers.h"
#include "config.h"

#define MPI_GET_UPPER_NEIGHBOUR(X) MPI_Cart_shift(grid, 0, 1, &X, &dist_junk);
#define MPI_GET_LOWER_NEIGHBOUR(X) MPI_Cart_shift(grid, 0, -1, &X, &dist_junk);
#define MPI_GET_LEFT_NEIGHBOUR(X) MPI_Cart_shift(grid, 1, 1, &X, &dist_junk);
#define MPI_GET_RIGHT_NEIGHBOUR(X) MPI_Cart_shift(grid, 1, -1, &X, &dist_junk);

#define SLEEPRANK std::this_thread::sleep_for(std::chrono::milliseconds(rank * 100))


void premultA(int row_size, int row_beginning, int K, int Ny, double *r0_l, double *x);

// Very specific parallel gradient descent just for this problem

void pcgd(  int rank, int size,     ///< MPI formulation
            double dx2, double dy2,    ///< Problem formulation
            int Ny, int K,          ///< A/b matrix formulation
            double* x, double* b);       