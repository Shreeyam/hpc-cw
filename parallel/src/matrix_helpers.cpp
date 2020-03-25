#include <iostream>
#include <fstream>
#include <iomanip>
#include "matrix_helpers.h"

using namespace std;

/* Conventions to be followed!
*  x is the horizontal axis, y is the vertical axis
*  i is the iterator for x, j is the iterator for y
*  Matrices are called in COLUMN MAJOR ormat
*  Thus, an index [i, j] -> [j + i * Ny]
*  Iterators to print row wise and then col wise will go
*  for(j) { for(i) { ... }}
*/

void printVec(int Ny, double* A)
{
    for (int j = 0; j < Ny; j++)
    {
        cout << A[j] << endl;
    }
}

void printVecH(int Ny, double* A)
{
    for (int j = 0; j < Ny; j++)
    {
        cout << A[j] << " ";
    }
    cout << endl;
}

void printVecH(int Ny, int* A)
{
    for (int j = 0; j < Ny; j++)
    {
        cout << A[j] << " ";
    }
    cout << endl;
}

void printMat(int Nx, int Ny, double* A)
{
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0; i < Nx; i++)
        {
            cout << A[j + i * Ny] << " ";
        }
        cout << endl;
    }
}

void dumpVec(int Ny, double* A)
{
    ofstream csv("dumpvec.csv", ios::out | ios::trunc);

    csv.precision(5);

    if (csv.good())
    {
        for (int j = 0; j < Ny; j++)
        {
            csv << A[j] << endl;
        }
    }

    csv.close();
}

void dumpMat(int Nx, int Ny, double* A)
{
    ofstream csv("dumpmatrix.csv", ios::out | ios::trunc);

    csv.precision(10);

    if (csv.good())
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                csv << A[j + i * Ny] << ", ";
            }
            csv << "," <<  endl;
        }
    }

    csv.close();
}

void zeroVec(int Ny, double* A)
{
    for (int j = 0; j < Ny; j++)
    {
        A[j] = 0;
    }
}

void zeroVec(int Ny, int* A)
{
    for (int j = 0; j < Ny; j++)
    {
        A[j] = 0;
    }
}

void zeroMat(int Nx, int Ny, double* A)
{
    zeroVec(Nx * Ny, A);
}

void copyVec(int Ny, double* A, double* B)
{
    for (int j = 0; j < Ny; j++)
    {
        B[j] = A[j];
    }
}

void copyMat(int Nx, int Ny, double* A, double*B)
{
    copyVec(Nx * Ny, A, B);
}

void addVec(int Ny, double* out, double* a, double* b, double beta)
{
    for(int i = 0; i < Ny; i++)
    {
        out[i] = a[i] + beta * b[i];
    }
}

double dot(int Ny, double *x, double *y)
{
    double sum = 0;

    for (int i = 0; i < Ny; i++)
    {
        sum += x[i] * y[i];
    }

    return sum;
}

double norm(int Ny, double* x)
{
    return dot(Ny, x, x);
}