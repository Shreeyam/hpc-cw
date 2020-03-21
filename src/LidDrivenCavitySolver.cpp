using namespace std;

#include "LidDrivenCavity.h"

#include <iostream>
#include "config.h"
#include "matrix_helpers.h"

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    cout << "Initialising solver..." << endl;
    // Configure the solver here...
    // ...
    solver->Initialise();

    // Run the solver
    solver->Integrate();

	return 0;
}
