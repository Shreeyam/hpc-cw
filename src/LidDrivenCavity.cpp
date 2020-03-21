#include "LidDrivenCavity.h"
#include "matrix_helpers.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    // Deallocaite pointers
    delete[] s;
    delete[] v;
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
    // I think that's it? Anything else...? Hmm...

    v = new double[Nx * Ny];
    s = new double[Nx * Ny];

    zeroMat(Nx, Ny, v);
    zeroMat(Nx, Ny, v);
}

void LidDrivenCavity::Integrate()
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

void LidDrivenCavity::updateBoundaries()
{
    // Update boundary conditions per side`
    // Top
    for (int i = 0; i < Nx; i++)
    {
        v[(Ny - 1) + i * Ny] = (s[(Ny - 1) + i * Ny] - s[(Ny - 2) + i * Ny]) * (2 / dy2) - ((2 * U) / (dy));
    }

    // Bottom
    for (int i = 0; i < Nx; i++)
    {
        v[(0) + i * Ny] = (s[(0) + i * Ny] - s[(1) + i * Ny]) * (2 / dy2);
    }

    // Left
    for (int j = 0; j < Ny; j++)
    {
        v[(j) + (0) * Ny] = (s[(j) + (0) * Ny] - s[(j) + (1) * Ny]) * (2/dx2);
    }

    // Right
    for (int j = 0; j < Ny; j++)
    {
        v[(j) + (0) * Ny] = (s[(j) + (Nx - 1) * Ny] - s[(j) + (Nx - 2) * Ny]) * (2/dx2);
    }
}

void LidDrivenCavity::updateInterior()
{

    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            v[i * Ny + j] = -(
                ((s[(i + 1) * Ny + j] - 2 * s[i * Ny + j] + s[(i - 1) * Ny + j]) / (dx2)) +
                ((s[i * Ny + (j + 1)] - 2 * s[i * Ny + j] + s[i * Ny + (j - 1)]) / (dy2)));
        }
    }
}
void LidDrivenCavity::newInterior()
{
}
void LidDrivenCavity::solvePoisson()
{
}