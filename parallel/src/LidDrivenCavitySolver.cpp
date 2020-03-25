#include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // for(int i = 0; i < argc;i++)
    // {

    // }

    // Initialise MPI.
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS)
    {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity *solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    solver->Initialise();

    // // Run the solver
    solver->Integrate();

	MPI_Finalize();
    return 0;
}
