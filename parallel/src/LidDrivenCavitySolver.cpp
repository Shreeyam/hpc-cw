#include <iostream>
#include <chrono>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Initialise MPI.
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS)
    {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity *solver = new LidDrivenCavity();

    int np;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
            Px = stoi(argv[10]);
            Py = stoi(argv[12]);
            dt = stod(argv[14]);
            T = stod(argv[16]);
            Re = stod(argv[18]);
        }
        catch (std::exception e)
        {
            if (rank == 0)
                cout << "Invalid arguments" << endl;

            return 0;
        }

        MPI_Comm_size(MPI_COMM_WORLD, &np);
        if (np != Px * Py)
        {
            if (rank == 0)
                cout << "np must equal Px * Py" << endl;

            return 0;
        }

        // If we make it here without any exceptions, print the arguments
        solver->SetDomainSize(Lx, Ly);
        solver->SetGridSize(Nx, Ny);
        solver->SetPartitions(Px, Py);
        solver->SetTimeStep(dt);
        solver->SetFinalTime(T);
        solver->SetReynoldsNumber(Re);

        if (rank == ROOTPROC)
        {
            cout << "Non-default parameters chosen:" << endl;
            cout << "Domain: " << Lx << " by " << Ly << " m " << endl;
            cout << "Grid cells: " << Nx << " by " << Ny << endl;
            cout << "Multicore partitions: " << Px << " by " << Py << endl;
            cout << "Timestep: " << dt << " s" << endl;
            cout << "Final time: " << T << " s" << endl;
            cout << "Re: " << Re << endl;
        }
    }

    // Configure the solver here...
    // ...
    if (rank == ROOTPROC)
        cout << "Initialising..." << endl;
    solver->Initialise();

    if (rank == ROOTPROC)
        cout << "Solving..." << endl;

    // Time it so the user knows
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Run the solver
    solver->Integrate();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    if (rank == ROOTPROC)
        std::cout << "Time elapsed:" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]" << std::endl;

    // LDC Destructor
    // Finalise MPI
    MPI_Finalize();
    return 0;
}
