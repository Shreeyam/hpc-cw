#include <thread>
#include <chrono>

// Central differences (in j and i)

#define CDi(x) (x[j + Ny * (i + 1)] - 2 * x[j + Ny * i] + x[j + Ny * (i - 1)])
#define CDj(x) (x[(j + 1) + Ny * (i)] - 2 * x[j + Ny * i] + x[(j - 1) + Ny * (i)])

// MPI Stuff for neighbours (assumes there's a junk variable defined)

#define MPI_GET_UPPER_NEIGHBOUR(X) MPI_Cart_shift(grid, 0, 1, &X, &junk);
#define MPI_GET_LOWER_NEIGHBOUR(X) MPI_Cart_shift(grid, 0, -1, &X, &junk);
#define MPI_GET_LEFT_NEIGHBOUR(X) MPI_Cart_shift(grid, 1, 1, &X, &junk);
#define MPI_GET_RIGHT_NEIGHBOUR(X) MPI_Cart_shift(grid, 1, -1, &X, &junk);

// Sleep 100 ms per rank
// Yeah it's hella stupid but at least it gets 
// cout to be sequential
#define SLEEPRANK std::this_thread::sleep_for(std::chrono::milliseconds(rank * 100))