#pragma once

#define DEBUG
// #define VERBOSE

// Default parameters
// Use const for better type-checking over #define...

// User parameters
const double LX = 1;
const double LY = 1;
const int NX = 12;
const int NY = 12;
const double PX = 2;
const double PY = 2;
const double TMAX = 1;
const double RE = 100;
const double DT = 1e-04;//0.01;

// Conjugate gradient descent parameters
const double CGD_TOL = 1E-18;
const double SIM_TOL = 1E-4;

// Not user parameters
// const doube