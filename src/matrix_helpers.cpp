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

    csv.precision(5);

    if (csv.good())
    {
        for (int j = 0; j < Ny; j++)
        {
            for (int i = 0; i < Nx; i++)
            {
                cout << A[j + i * Ny] << " ";
            }
            cout << "," <<  endl;
        }
    }

    csv.close();
}
