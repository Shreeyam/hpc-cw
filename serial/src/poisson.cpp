#include "poisson.h"

Poisson::Poisson()
{

}

void Poisson::constructA(int Nx, int Ny, double dx2, double dy2)
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

        if ((i + 1) % (Ny - 2) != 0)
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

Poisson::~Poisson()
{
    delete[] P;
    delete[] A;
}

void Poisson::solvePoisson(int Nx, int Ny, double dx2, double dy2, double *v, double *s)
{
    // dgbsv!
    // KL = KU = Ny
    // For parallelisation, look at pdgbsv (Example 26.2) in the notes
    // Solve A\psi = \omega = As = v
    // A is overwritten by LU factorisation, so needs to be constructed again

    A = new double[((Nx - 2) * (Ny - 2)) * (3 * (Ny - 2) + 1)];
    P = new int[(Nx-2) * (Ny-2)];

    zeroVec((Nx-2) * (Ny-2), P);

    constructA(Nx, Ny, dx2, dy2);

    double *b = new double[(Nx - 2) * (Ny - 2)];
    // Constructed A, now have to construct a working b matrix
    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            b[(j - 1) + (i - 1) * (Ny - 2)] = v[j + i * Ny];
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
        (Nx - 2) * (Ny - 2), ///< LDB
        info                 ///< INFO
    );


    // Transfer answer into s, leave boundaries alone
    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            s[j + i * Ny] = b[(j - 1) + (i - 1) * (Ny - 2)];
        }
    }

    delete[] b;
    delete[] A;
    delete[] P;

    // cout << "New s:" << endl;
    // printMat(Nx, Ny, s);

#ifdef VERBOSE
    cout << "info:" << info << endl;
#endif
}
