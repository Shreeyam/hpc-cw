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

    v = new double[Nx * Ny];
    s = new double[Nx * Ny];
    b = new double[(Nx - 2) * (Ny - 2)];

    A = new double[((Nx - 2) * (Ny - 2)) * (3 * (Ny - 2) + 1)];
    P = new int[Nx * Ny];

#ifdef VERBOSE
    cout << "Allocated memory..." << endl;
#endif

    zeroMat(Nx, Ny, v);
    zeroMat(Nx, Ny, s);

    zeroVec(Nx * Ny, P);

    dt = DT_MAX;

#ifdef VERBOSE
    cout << "Zeroed matrices..." << endl;
    cout << "dx: " << dx << endl;
    cout << "dy: " << dy << endl;
// cout << "dt: " << dt << endl;
#endif
}

void LidDrivenCavity::Integrate()
{
#ifdef VERBOSE
    cout << "Solving..." << endl;
#endif

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

        cout << i << endl;

#ifdef VERBOSE
        cout << "Vorticity:" << endl;
        printMat(Nx, Ny, v);
        cout << "Stream:" << endl;
        printMat(Nx, Ny, s);
#endif
    }

    dumpMat(Nx, Ny, s);
}

void LidDrivenCavity::updateBoundaries()
{
    // Update boundary conditions per side`
    // Top
    // cout << "Pre-anything";
    // printMat(Nx, Ny, v);

    for (int i = 0; i < Nx; i++)
    {
        v[(Ny - 1) + i * Ny] = (s[(Ny - 1) + i * Ny] - s[(Ny - 2) + i * Ny]) * (2 / dy2) - ((2 * U) / (dy));
    }

    // cout << "Post-top BC";
    // printMat(Nx, Ny, v);

    // Bottom
    for (int i = 0; i < Nx; i++)
    {
        v[(0) + i * Ny] = (s[(0) + i * Ny] - s[(1) + i * Ny]) * (2 / dy2);
    }

    // cout << "Post-bottom BC";
    // printMat(Nx, Ny, v);

    // Left
    for (int j = 0; j < Ny; j++)
    {
        v[(j) + (0) * Ny] = (s[(j) + (0) * Ny] - s[(j) + (1) * Ny]) * (2 / dx2);
    }

    // cout << "Post-left BC";
    // printMat(Nx, Ny, v);

    // Right
    for (int j = 0; j < Ny; j++)
    {
        v[(j) + (Nx - 1) * Ny] = (s[(j) + (Nx - 1) * Ny] - s[(j) + (Nx - 2) * Ny]) * (2 / dx2);
    }

    // cout << "Post-right BC";
    // printMat(Nx, Ny, v);
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
            v[i * Ny + j] += dt * (
                -1 * (((s[i * Ny + (j + 1)] - s[i * Ny + (j - 1)]) / (2 * dy)) * ((v[(i + 1) * Ny + (j)] - v[(i - 1) * Ny + (j)]) / (2 * dx)))
                   + (((s[(i + 1) * Ny + j] - s[(i - 1) * Ny + j]) / (2 * dx)) * ((v[(i)*Ny + (j + 1)] - v[(i)*Ny + (j - 1)]) / (2 * dy))) + 
                   ((1 / Re) * 
                   (((v[(i + 1) * Ny + (j)] - 2 * v[(i)*Ny + (j)] + v[(i - 1) * Ny + (j)]) / (dx2)) + ((v[(i)*Ny + (j + 1)] - 2 * v[(i)*Ny + (j)] + v[(i)*Ny + (j - 1)]) / (dy2)))));
        }
    }
}

void LidDrivenCavity::solvePoisson()
{
    // dgbsv!
    // KL = KU = Ny
    // For parallelisation, look at pdgbsv (Example 26.2) in the notes
    // Solve A\psi = \omega = As = v
    // A is overwritten by LU factorisation, so needs to be constructed again

    constructA();

    // Constructed A, now have to construct a working b matrix
    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            b[(j-1) + (i-1) * (Ny-2)] = v[j + i * Ny];
        }
    }

    /// Copy v into s so we don't need to allocate any extra working memory
    // copyMat(Nx, Ny, v, s);

    int info = 0;

    // God these Doxygen comments are awful.
    // I thought the point of comments was to be easy,
    // not to write 4 characters out each time I want
    // to comment a single line...

    F77NAME(dgbsv)
    (
        (Nx - 2) * (Ny - 2), ///< N
        (Ny - 2),            ///< KL
        (Ny - 2),            ///< KU
        1,                   ///< NRHS
        A,                   ///< AB
        3 * (Ny - 2) + 1,    ///< LDAB
        P,                   ///< IPIV
        b,                   ///< B
        (Nx-2) * (Ny-2),             ///< LDB
        info                 ///< INFO
    );

    // Transfer answer into s, leave boundaries alone
    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            s[j + i * Ny] = b[(j-1) + (i-1) * (Ny-2)];
        }
    }

    // cout << "New s:" << endl;
    // printMat(Nx, Ny, s);

    #ifdef VERBOSE
    cout << "info:" << info << endl;
    #endif
}

void LidDrivenCavity::constructA()
{
    // Split each term on a new line or I'm never gonna be able to read this
    zeroMat((Nx - 2) * (Ny - 2), 3 * (Ny - 2) + 1, A);

    for (int i = 0; i < (Nx - 2) * (Ny - 2); i++)
    {
        // Offsets of +1 for diagonal and -1 for 
        // matrix orientation cancel eachother out
        if ((i) % (Ny - 2) != 0)
        {
            // j = Ny-1 (first superdiagonal)
            A[i * (3 * (Ny - 2) + 1) + ((Ny - 2) - 1 + (Ny - 2))] = -1 / dy2;
        }

        if ((i+1) % (Ny - 2) != 0)
        {
            A[i * (3 * (Ny - 2) + 1) + (2 * (Ny - 2) + 1)] = -1 / dy2;
        }
        // j = 0 (highest superdiagonal)
        A[i * (3 * (Ny - 2) + 1) + ((Ny - 2))] = -1 / dx2;
        // j = Ny (leading diagonal)
        A[i * (3 * (Ny - 2) + 1) + (2 * (Ny - 2))] = 2 / dy2 + 2 / dx2;

        A[i * (3 * (Ny - 2) + 1) + (3 * (Ny - 2))] = -1 / dx2;
    }
}