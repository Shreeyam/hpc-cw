#include "LidDrivenCavity.h"
#include "matrix_helpers.h"
#include "poisson.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    // Deallocaite pointers
    delete[] s;
    delete[] v;
    delete[] b;
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

    v = new double[Nx * Ny];
    s = new double[Nx * Ny];
    b = new double[(Nx - 2) * (Ny - 2)];


#ifdef VERBOSE
    cout << "Allocated memory..." << endl;
#endif

    zeroMat(Nx, Ny, v);
    zeroMat(Nx, Ny, s);

    poisson = new Poisson();

#ifdef VERBOSE
    cout << "Zeroed matrices..." << endl;
    cout << "dx: " << dx << endl;
    cout << "dy: " << dy << endl;
// cout << "dt: " << dt << endl;
#endif
}

void LidDrivenCavity::Integrate()
{
    for (int i = 0; i < (T / dt) + 1; i++)
    {
        // Steps!

        // 1. Update v at boundaries
        updateBoundaries();
        // 2. Update v in interior
        updateInterior();
        // 3. Compute vorticity at t + dt in the interior
        newInterior();
        // 4. Solve poisson equation to compute s at time t + dt
        poisson->solvePoisson(Nx, Ny, dx2, dy2, v, s);

        cout << i << endl;

#ifdef VERBOSE
        cout << "Vorticity:" << endl;
        printMat(Nx, Ny, v);
        cout << "Stream:" << endl;
        printMat(Nx, Ny, s);
#endif
    }

    // Output to file
    dumpMat(Nx, Ny, v, (char *)"vorticity.csv");
    dumpMat(Nx, Ny, s, (char *)"streamfunc.csv");
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
        v[(j) + (0) * Ny] = (s[(j) + (0) * Ny] - s[(j) + (1) * Ny]) * (2 / dx2);
    }

    // Right
    for (int j = 0; j < Ny; j++)
    {
        v[(j) + (Nx - 1) * Ny] = (s[(j) + (Nx - 1) * Ny] - s[(j) + (Nx - 2) * Ny]) * (2 / dx2);
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
    // Split each term on a new line or I'm never gonna be able to read this
    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            double delta = (-1 * (((s[i * Ny + (j + 1)] - s[i * Ny + (j - 1)]) / (2 * dy)) * ((v[(i + 1) * Ny + (j)] - v[(i - 1) * Ny + (j)]) / (2 * dx))) + (((s[(i + 1) * Ny + j] - s[(i - 1) * Ny + j]) / (2 * dx)) * ((v[(i)*Ny + (j + 1)] - v[(i)*Ny + (j - 1)]) / (2 * dy))) +
                            ((1 / Re) *
                             (((v[(i + 1) * Ny + (j)] - 2 * v[(i)*Ny + (j)] + v[(i - 1) * Ny + (j)]) / (dx2)) + ((v[(i)*Ny + (j + 1)] - 2 * v[(i)*Ny + (j)] + v[(i)*Ny + (j - 1)]) / (dy2)))));
            v[i * Ny + j] += dt * delta;
        }
    }
}
