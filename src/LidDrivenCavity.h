#pragma once

#include <string>
#include <math.h>
#include "config.h" 

using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);

    void Initialise();
    void Integrate();

    // Add any other public functions

private:
    double* v = nullptr;    ///< Vorticity
    double* s = nullptr;    ///< Stream function

    double* A = nullptr;    ///< Poisson equation matrix
    // DPBTRS
    // User parameters (initialise to defaults)
    double dt = DT;
    double T = T;
    int    Nx = NX;
    int    Ny = NY;
    double Lx = LX;
    double Ly = LY;
    double Re = RE;

    // Non user-modifiable parameters
    const double U = 1.0;

    // Derived parameters
    const double dx = Lx/(Nx - 1);
    const double dy = Ly/(Nx - 1);
    const double DT_MAX = (Re * dx * dy) / 4;

    // Precomputed powers
    const double dx2 = 1;//pow(dx, 2);  
    const double dy2 = 1;//pow(dy, 2);  

    void updateBoundaries();
    void updateInterior();
    void newInterior();
    void solvePoisson();
    void constructA();
    
};

