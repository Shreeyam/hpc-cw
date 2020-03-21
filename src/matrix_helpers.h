#include <iostream>

/*
* @brief Dumps a vector to stdout
* @param Ny The length of the vector in total element count
* @param A The pointer to the vector itself
*/
void printVec(int Ny, double* A);

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

/*
* @brief Zeros all the elements of a matrix 
* @param Nx is the length of the vector in elements in the x axis
* @param Ny The length of the vector in elements in the y axis
* @param A The pointer to the vector itself
*/
void zeroMat(int Nx, int Ny, double* A);