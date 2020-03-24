#include <iostream>

/*
* @brief Dumps a vector to stdout
* @param Ny The length of the vector in total element count
* @param A The pointer to the vector itself
*/
void printVec(int Ny, double* A);
void printVecH(int Ny, double* A);

/*
* @brief Dumps a (col-major!) matrix to stdout 
* @param Nx is the length of the vector in elements in the x axis
* @param Ny The length of the vector in elements in the y axis
* @param A The pointer to the vector itself
*/
void printMat(int Nx, int Ny, double* A);

/*
* @brief Dumps a vector to a csv file
* @param Ny The length of the vector in total element count
* @param A The pointer to the vector itself
*/
void dumpVec(int Ny, double* A);

/*
* @brief Dumps a (col-major!) matrix to a csv file 
* @param Nx is the length of the vector in elements in the x axis
* @param Ny The length of the vector in elements in the y axis
* @param A The pointer to the vector itself
*/
void dumpMat(int Nx, int Ny, double* A);

/*
* @brief Zeros all the elements of a vector
* @param Ny The length of the vector in total element count
* @param A The pointer to the vector itself
*/
void zeroVec(int Ny, double* A);
void zeroVec(int Ny, int* A);

/*
* @brief Zeros all the elements of a matrix 
* @param Nx is the length of the vector in elements in the x axis
* @param Ny The length of the vector in elements in the y axis
* @param A The pointer to the vector itself
*/
void zeroMat(int Nx, int Ny, double* A);

// Copies A into B
void copyVec(int Ny, double* A, double*B);
void copyMat(int Nx, int Ny, double* A, double*B);

// out <- a + \beta b
void addVec(int Ny, double* out, double* a, double* b, double beta);

double dot(int Ny, double* x, double* y);
double norm(int Ny, double* x);