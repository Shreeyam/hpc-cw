#include "LidDrivenCavity.h"
#include "matrix_helpers.h"

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(dgbsv)(
        const int &n,
        const int &kl,
        const int &ku,
        const int &nrhs,
        const double *AB,
        const int &ldab,
        int *ipiv,
        double *B,
        const int &ldb,
        int &info);
}

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    // Deallocaite pointers
    delete[] s;
    delete[] v;
    delete[] b;

    delete[] A;
    delete[] P;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx = nx;
    Ny = ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re = re;
}

void LidDrivenCavity::Initialise()
{
    // TODO:
    // Assign the memory for the vorticity, stream function
    // Initialise that memory to all zeros
    // Precompute an A matrix (for the Poisson solver)
    // I think that's it? Anything else...? Hmm...

    // Steps I need to do here:

    
    v = new double[Nx * Ny];
    s = new double[Nx * Ny];


}

void LidDrivenCavity::Integrate()
{

    for (int i = 0; i < 1000; i++)
    {
        // Steps!

        // 1. Update v at boundaries
        updateBoundaries();
        // 2. Update v in interior
        updateInterior();
        // 3. Compute vorticity at t + dt in the interior
        newInterior();
        // 4. Solve poisson equation to compute s at time t + dt
        solvePoisson();

    }
}

void LidDrivenCavity::updateBoundaries()
{

}

void LidDrivenCavity::updateInterior()
{

}
void LidDrivenCavity::newInterior()
{

}

void LidDrivenCavity::solvePoisson()
{

}
