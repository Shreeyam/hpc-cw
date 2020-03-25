#pragma once

#include <iostream>
#include <math.h>
#include <string>
#include "matrix_helpers.h"
#include "parallel.h"
#include "macros.h"
#include "config.h"

using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetPartitions(int px, int py);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);

    void Initialise();
    void Integrate();

    // Add any other public functions

private:
    // Global variables, not allocated except for output/exchange on root!
    double *v = nullptr; ///< Vorticity
    double *s = nullptr; ///< Stream function

    // Local variables, with guard halos
    double *v_l = nullptr;
    double *s_l = nullptr;

    // Serialised regions
    double *v_s = nullptr;
    double *s_s = nullptr;

    /// User parameters (initialise to defaults fron config.h)
    double dt = DT;
    double T = TMAX;
    int Nx = NX;
    int Ny = NY;
    double Lx = LX;
    double Ly = LY;
    double Re = RE;

    // Parallel computing
    int rank;
    int sizes[2];
    int periods[2] = {0, 0};
    int reorder = 0;
    int coords[2];
    MPI_Comm grid;

    // Neighbour ranks
    int neighbour_up;
    int neighbour_down;
    int neighbour_left;
    int neighbour_right;

    // Parallel partitions
    int Px = PX;
    int Py = PY;

    // Average block sizes + remainders
    int Nyb;
    int Nxb;

    int rem_y;
    int rem_x;

    // Local block sizes (of THIS process)
    int Nyl;
    int Nxl;

    // Skip sizes (instead of +2 each time...)
    int Nylh;

    // Bounds for interior calculations
    int ibounds_j0;
    int ibounds_j1;
    int ibounds_i0;
    int ibounds_i1;

    // Total number of processes
    int P = Px * Py;

    // Data buffers for receiving
    double *txtopBuffer;
    double *txbottomBuffer;
    double *txleftBuffer;
    double *txrightBuffer;

    double *rxtopBuffer;
    double *rxbottomBuffer;
    double *rxleftBuffer;
    double *rxrightBuffer;

    double *txrxBuffer;

    bool converged = false;

    int *txrxcounts;
    int *displs;

    // Discardable parameters
    int junk;

    /// Non user-modifiable parameters
    const double U = 1.0;

    /// Derived parameters
    const double dx = Lx / (Nx - 1);
    const double dy = Ly / (Nx - 1);
    const double DT_MAX = (Re * dx * dy) / 4;

    /// Precomputed powers
    const double dx2 = pow(dx, 2);
    const double dy2 = pow(dy, 2);

    void updateHalo();
    void updatevHalo();
    void updatesHalo();
    void txrx();

    /// Solver steps
    void updateBoundaries();
    void updateInterior();
    void newInterior();
    void solvePoisson();

    /// Unpacks piecewise matrices to the full thing
    /// (only on process 0 for output)
    void filltxbuf(double* source);
    void emptytxbuf(double* dest);
    void deserialize(double *v_s, double *v);
    void serialize(double *v_s, double *v);
};
