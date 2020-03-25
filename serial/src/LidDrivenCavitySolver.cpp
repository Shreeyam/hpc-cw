using namespace std;

#include "LidDrivenCavity.h"

#include <iostream>
#include <chrono>
#include "config.h"
#include "matrix_helpers.h"

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity *solver = new LidDrivenCavity();

    // Skip argument 1
    // If we're not using the defaults...
    if (argc > 1)
    {
        // Some sanity checks...
        double Lx, Ly, dt, T, Re;
        int Nx, Ny, Px, Py;

        try
        {
            // Try and do some conversions
            Lx = stod(argv[2]);
            Ly = stod(argv[4]);
            Nx = stoi(argv[6]);
            Ny = stoi(argv[8]);
            dt = stod(argv[10]);
            T = stod(argv[12]);
            Re = stod(argv[14]);
        }
        catch (std::exception e)
        {
            cout << "Invalid arguments" << endl;
        }

        // If we make it here without any exceptions, print the arguments
        solver->SetDomainSize(Lx, Ly);
        solver->SetGridSize(Nx, Ny);
        solver->SetTimeStep(dt);
        solver->SetFinalTime(T);
        solver->SetReynoldsNumber(Re);

        cout << "Non-default parameters chosen:" << endl;
        cout << "Domain: " << Lx << " by " << Ly << " m " << endl;
        cout << "Grid cells: " << Nx << " by " << Ny << endl;
        cout << "Timestep: " << dt << " s" << endl;
        cout << "Final time: " << T << " s" << endl;
        cout << "Re: " << Re << endl;
    }

    // Configure the solver here...
    // ...
    solver->Initialise();

    // Time it so the user knows
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Run the solver
    solver->Integrate();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time elapsed:" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]" << std::endl;

    // LDC Destructor
    // Finalise MPI
    return 0;
}
