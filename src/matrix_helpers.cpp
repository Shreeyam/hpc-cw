#include <iostream>
#include <fstream>
#include <iomanip>
#include "matrix_helpers.h"

using namespace std;

// Conventions to be followed!
// x is the horizontal axis, y is the vertical axis
// i is the iterator for x, j is the iterator for y
// Matrices are called in COLUMN MAJOR ormat
// Thus, an index [i, j] -> [j + i * Ny]
// Iterators to print row wise and then col wise will go
// for(j) { for(i) { ... }}

/*
* Dumps a vector to stdout
* Ny is the length of the vector in total element count
* A is the pointer to the vector itself
*/
void printVec(int Ny, double* A)
{
    for (int j = 0; j < Ny; j++)
    {
        cout << A[j] << endl;
    }
}

/*
* Dumps a (col-major!) matrix to stdout 
* Nx is the length of the vector in elements in the x axis
* Ny is the length of the vector in elements in the y axis
* A is the pointer to the vector itself
*/
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

/*
* Dumps a vector to a csv file
* Ny is the length of the vector in total element count
* A is the pointer to the vector itself
*/
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

/*
* Dumps a (col-major!) matrix to a csv file 
* Nx is the length of the vector in elements in the x axis
* Ny is the length of the vector in elements in the y axis
* A is the pointer to the vector itself
*/
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
