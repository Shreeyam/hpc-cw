#pragma once

#include <string>
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
    double* v = nullptr;    // Vorticity
    double* s = nullptr;    // Stream function

    // User parameters
    double dt = DT;
    double T = T;
    int    Nx = NX;
    int    Ny = NY;
    double Lx = LX;
    double Ly = LY;
    double Re = RE;

    // Non user-modifiable parameters

    const double DX = Lx/(Nx - 1);
    const double DY = Ly/(Nx - 1);
    const double DT_MAX = (Re * DX * DY) / 4;
    
};

