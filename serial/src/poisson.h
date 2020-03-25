#pragma once
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

using namespace std;

class Poisson
{
public:
    Poisson();
    ~Poisson();
    void solvePoisson(int Nx, int Ny, double dx2, double dy2, double *v, double *s);

private:
    void constructA(int Nx, int Ny, double dx2, double dy2);
    int *P = nullptr; ///< Pivot vector
    double *A = nullptr;
};