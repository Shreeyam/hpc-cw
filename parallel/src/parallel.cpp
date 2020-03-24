#include "parallel.h"

int dist_junk;

void premultA(int row_size, int row_beginning, int K, int Ny, double *r0_l, double *x, double *b_l, double alpha)
{
    for (int i = 0; i < row_size; i++)
    {
        // Actual row on the matrix
        int j = i + row_beginning;
        r0_l[i] = 64 * x[j];

        if (j + 1 < Ny && (j + 1) % (K) != 0)
        {
            r0_l[i] += -16 * x[j + 1];
        }

        if (j + K < Ny)
        {
            r0_l[i] += -16 * x[j + K];
        }

        if (j >= K)
        {
            r0_l[i] += -16 * x[j - K];
        }

        if (j >= 1 && (j) % (K) != 0)
        {
            r0_l[i] += -16 * x[j - 1];
        }

        if (b_l != nullptr)
        {
            r0_l[i] = b_l[i] + alpha * r0_l[i];
        }
    }
}

void pcgd(int rank, int size,     ///< MPI formulation
          double dx2, double dy2, ///< Problem formulation
          int Ny, int K,          ///< A/b matrix formulation
          double *x, double *b)
{
    // Helper variable when referring to size
    // in MPI contexts and P in matrix contexts
    const int P = size;

    // Average length to split (local)
    const int Nyl = (Ny / P);
    // Remainder on y that don't divide
    const int rem = (Ny % P);

    int Nbl = (rank != size - 1) ? Nyl : Nyl + rem;

    // MPI Setup for data sharing
    int *sendcounts_y = new int[size];
    int *displs_y = new int[size];

    // Divide the data
    for (int i = 0; i < size; i++)
    {
        if (i != size - 1)
        {
            sendcounts_y[i] = Nyl;
        }
        else
        {
            sendcounts_y[i] = Nyl + rem;
        }

        displs_y[i] = Nyl * i;
    }

    // Create localised storage of b matrices
    double *b_l = new double[Nbl];
    MPI_Scatterv(b, sendcounts_y, displs_y, MPI_DOUBLE, b_l, Nbl, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Now we have a local b, need to construct a local x... and finish with the multiply... Just need to calculate b - Ax0

    // Figure out local row beginning... (0 indexed)
    // 0s every K entires on the diagonal with an offset

    // On each process...
    int row_beginning = Nyl * rank;
    // What we're all computing...
    int row_size = Nbl;

    double *r0 = new double[Ny];
    double *r1 = new double[Ny];
    double *p0 = new double[Ny];
    double *r0_l = new double[Nbl];
    double *r1_l = new double[Nbl];
    double *p0_l = new double[Nbl];

    double *Apk_l = new double[Nbl];

    double pApk_l;
    double pApk;
    double r0norm_l;
    double r0norm;

    int k = 0;
    // alpha in cgd
    double ak = 0;
    double bk = 0;

    zeroVec(Nbl, r0_l);

    // Start of conjugate gradient descent...

    double error;

    // Matrix multiplication and subtraction for r0_l and x
    // Lucky fucking guess at the indices
    premultA(row_size, row_beginning, K, Ny, r0_l, x, b_l, -1.0);

    copyVec(Nbl, r0_l, p0_l);

    MPI_Allgatherv(r0_l, Nbl, MPI_DOUBLE, r0, sendcounts_y, displs_y, MPI_DOUBLE, MPI_COMM_WORLD);

    // Collected to p0
    copyVec(Ny, r0, p0);

    // Return back to parallel operations...

    do
    {
        premultA(row_size, row_beginning, K, Ny, Apk_l, p0, nullptr, -1.0);

        // Don't forget to redistribute locals!

        pApk_l = dot(Nbl, p0_l, Apk_l);
        r0norm_l = norm(Nbl, r0_l);

        MPI_Allreduce(&pApk_l, &pApk, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&r0norm_l, &r0norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Got the dot products Now
        ak = r0norm / pApk;

        addVec(Ny, x, x, p0, ak);
        addVec(Nbl, r1_l, r0_l, Apk_l, -ak);

        MPI_Allgatherv(r1_l, Nbl, MPI_DOUBLE, r1, sendcounts_y, displs_y, MPI_DOUBLE, MPI_COMM_WORLD);

        // Calculate error
        error = sqrt(norm(Ny, r1));

        // More dot products... again for denominator and
        double rk1dot_l = norm(Nbl, r1_l);
        double rk0dot_l = norm(Nbl, r0_l);

        // Store from all other processes here
        double rk1dot;
        double rk0dot;

        // Reduce to a beta term
        MPI_Allreduce(&rk1dot_l, &rk1dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&rk0dot_l, &rk0dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        bk = rk1dot / rk0dot;

        addVec(Nbl, p0_l, r1_l, p0_l, bk);

        MPI_Allgatherv(p0_l, Nbl, MPI_DOUBLE, p0, sendcounts_y, displs_y, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(r1_l, Nbl, MPI_DOUBLE, r0, sendcounts_y, displs_y, MPI_DOUBLE, MPI_COMM_WORLD);

        copyVec(Nbl, r1_l, r0_l);

        k++;

    } while (error > CGD_TOL);

    if (rank == 0)
    {
        printVecH(Ny, x);
    }

    delete[] r0;
    delete[] r1;
    delete[] p0;
    delete[] r0_l;
    delete[] r1_l;
    delete[] Apk_l;
    delete[] sendcounts_y;
    delete[] displs_y;
    delete[] b_l;
}