#include "LidDrivenCavity.h"

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
    // Deallocaite pointers
    delete[] s;
    delete[] v;
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{
    Lx = xlen;
    Ly = ylen;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
    Nx = nx;
    Ny = ny;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
    dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
    T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
    Re = re;
}

void LidDrivenCavity::Initialise()
{
    // TODO:
    // Assign the memory for the vorticity, stream function
    // Initialise that memory to all zeros
    // Precompute an A matrix (for the Poisson solver)
    // I think that's it? Anything else...? Hmm...

    // Steps I need to do here:
    // Get in the mindset, I am a process
    // I am allocated a block with a bunch of guard cells to do whatever with
    // First order of business: Where am I? on the cartesian grid
    // How many elements do I have top and bottom of me?

    v = new double[Nx * Ny];
    s = new double[Nx * Ny];

    // Serial buffer for process 0
    v_s = new double[Nx * Ny];
    s_s = new double[Nx * Ny];

    // What's the total size?
    sizes[0] = Py;
    sizes[1] = Px;

    MPI_Cart_create(MPI_COMM_WORLD, 2, sizes, periods, reorder, &grid);

    // Get rank and co-ordinates
    MPI_Comm_rank(grid, &rank);
    MPI_Cart_coords(grid, rank, 2, coords);

    // Find my neighbours
    MPI_GET_UPPER_NEIGHBOUR(neighbour_up);
    MPI_GET_LOWER_NEIGHBOUR(neighbour_down);
    MPI_GET_LEFT_NEIGHBOUR(neighbour_left);
    MPI_GET_RIGHT_NEIGHBOUR(neighbour_right);

    // What's the normal block size?
    Nyb = Ny / Py;
    Nxb = Nx / Py;

    // What's left over?
    rem_y = Ny % Py;
    rem_x = Nx % Px;

    // What's my block size?
    // Am I on the right or lower boundaries?
    Nyl = (coords[0] == Py - 1) ? Nyb + rem_y : Nyb;
    Nxl = (coords[1] == Px - 1) ? Nxb + rem_x : Nxb;

    // Local skip
    Nylh = Nyl + 2;

    // Allocate my local memory... (with halo)
    v_l = new double[(Nyl + 2) * (Nxl + 2)];
    s_l = new double[(Nyl + 2) * (Nxl + 2)];

    // ...and zero it out as per the initial conditions
    zeroMat(Nxl + 2, Nyl + 2, v_l);
    zeroMat(Nxl + 2, Nyl + 2, s_l);

    // Allocate some space to each buffer
    txtopBuffer = new double[Nxl];
    txbottomBuffer = new double[Nxl];
    txleftBuffer = new double[Nyl];
    txrightBuffer = new double[Nyl];

    rxtopBuffer = new double[Nxl];
    rxbottomBuffer = new double[Nxl];
    rxleftBuffer = new double[Nyl];
    rxrightBuffer = new double[Nyl];

    // Serial buffers for step 4
    txrxBuffer = new double[Nxl * Nyl];

    // Calculate bounds
    // If on top boundary => j0 = 1, else 0
    // If on bottom boundary => j1 = Nyl - 1, else Nxl
    // If on left boundary => i0 = 1, else 0
    // If on right boundary => i1 = Nxl - 1, else Nxl
    ibounds_j0 = (coords[0] == 0) ? 1 : 0;
    ibounds_j1 = (coords[0] == Py - 1) ? Nyl - 1 : Nyl;
    ibounds_i0 = (coords[1] == 0) ? 1 : 0;
    ibounds_i1 = (coords[1] == Px - 1) ? Nxl - 1 : Nxl;

    // Fill up send counts and displacements

    txrxcounts = new int[P];
    displs = new int[P];

    zeroVec(P, displs);

    for (int j = 0; j < Py; j++)
    {
        for (int i = 0; i < Px; i++)
        {
            if (i == Px - 1 && j == Py - 1)
            {

                txrxcounts[j + Py * i] = (Nyb + rem_y) * (Nxb + rem_x);
            }
            else if (i == Px - 1)
            {
                txrxcounts[j + Py * i] = (Nyb) * (Nxb + rem_x);
            }
            else if (j == Py - 1)
            {
                txrxcounts[j + Py * i] = (Nyb + rem_y) * (Nxb);
            }
            else
            {
                txrxcounts[j + Py * i] = (Nyb) * (Nxb);
            }
        }
    }

    for (int i = 1; i < P; i++)
    {
        displs[i] = displs[i - 1] + txrxcounts[i - 1];
    }
}

void LidDrivenCavity::Integrate()
{
    for (int i = 0; i < T / dt + 1; i++)
    {
        // Steps!
        updateHalo();
        // 1. Update v at boundaries

        updateBoundaries();

        // 2. Update v in interior
        updateInterior();

        // 3. Compute vorticity at t + dt in the interior
        newInterior();

        // // 4. Solve poisson equation to compute s at time t + dt
        solvePoisson();
    }

    if (rank == ROOTPROC)
    {
        printMat(Nx, Ny, s);
    }
}

void LidDrivenCavity::updateBoundaries()
{
    // Am I a boundary process? Do I need to care?
    // Top boundary
    if (coords[0] == Py - 1)
    {
        for (int i = 0; i < Nxl; i++)
        {
            // New coordinates: k = i+1, l = j + 1
            // +1/+1 offset to account for halo cells
            v_l[(Nyl - 1 + 1) + (i + 1) * Nylh] = (s_l[(Nyl - 1 + 1) + (i + 1) * Nylh] - s_l[(Nyl - 2 + 1) + (i + 1) * Nylh]) * (2 / dy2) - ((2 * U) / (dy));
        }
    }

    // Bottom boundary
    if (coords[0] == 0)
    {
        for (int i = 0; i < Nxl; i++)
        {
            v_l[(0 + 1) + (i + 1) * Nylh] = (s_l[(0 + 1) + (i + 1) * Nylh] - s_l[(1 + 1) + (i + 1) * Nylh]) * (2 / dy2);
        }
    }

    // Left boundary
    if (coords[1] == 0)
    {
        for (int j = 0; j < Nyl; j++)
        {
            v_l[(j + 1) + (0 + 1) * Nylh] = (s_l[(j + 1) + (0 + 1) * Nylh] - s_l[(j + 1) + (1 + 1) * Nylh]) * (2 / dx2);
        }
    }

    // Right boundary
    if (coords[1] == Px - 1)
    {
        for (int j = 0; j < Nyl; j++)
        {
            v_l[(j + 1) + (Nxl - 1 + 1) * Nylh] = (s_l[(j + 1) + (Nxl - 1 + 1) * Nylh] - s_l[(j + 1) + (Nxl - 2 + 1) * Nylh]) * (2 / dx2);
        }
    }
}

void LidDrivenCavity::updateInterior()
{
    // Calculating interiors... Need to be careful about limits

    for (int j = ibounds_j0; j < ibounds_j1; j++)
    {
        for (int i = ibounds_i0; i < ibounds_i1; i++)
        {
            v[i * Ny + j] = -(
                ((s[(i + 1 + 1) * Nylh + (j + 1)] - 2 * s[(i + 1) * Nylh + (j + 1)] + s[(i - 1 + 1) * Nylh + (j + 1)]) / (dx2)) +
                ((s[(i + 1) * Nylh + (j + 1 + 1)] - 2 * s[(i + 1) * Nylh + (j + 1)] + s[(i + 1) * Nylh + (j - 1 + 1)]) / (dy2)));
        }
    }
}
void LidDrivenCavity::newInterior()
{
    // Split each term on a new line or I'm never gonna be able to read this
    // This is a mess because I did a find and replace on the serial code
    // But you think I'd touch it if it works? No way
    double deltamax = 0;

    for (int j = ibounds_j0; j < ibounds_j1; j++)
    {
        for (int i = ibounds_i0; i < ibounds_i1; i++)
        {
            double delta = (-1 * (((s_l[(i + 1) * Nylh + ((j + 1) + 1)] - s_l[(i + 1) * Nylh + ((j + 1) - 1)]) / (2 * dy)) * ((v_l[((i + 1) + 1) * Nylh + ((j + 1))] - v_l[((i + 1) - 1) * Nylh + ((j + 1))]) / (2 * dx))) + (((s_l[((i + 1) + 1) * Nylh + (j + 1)] - s_l[((i + 1) - 1) * Nylh + (j + 1)]) / (2 * dx)) * ((v_l[((i + 1)) * Nylh + ((j + 1) + 1)] - v_l[((i + 1)) * Nylh + ((j + 1) - 1)]) / (2 * dy))) +
                            ((1 / Re) *
                             (((v_l[((i + 1) + 1) * Nylh + ((j + 1))] - 2 * v_l[((i + 1)) * Nylh + ((j + 1))] + v_l[((i + 1) - 1) * Nylh + ((j + 1))]) / (dx2)) + ((v_l[((i + 1)) * Nylh + ((j + 1) + 1)] - 2 * v_l[((i + 1)) * Nylh + ((j + 1))] + v_l[((i + 1)) * Nylh + ((j + 1) - 1)]) / (dy2)))));

            v_l[(i + 1) * Nylh + (j + 1)] += dt * delta;

            deltamax = max(deltamax, abs(delta));
        }
    }

    if (deltamax < SIM_TOL)
    {
        converged = true;
    }

    // cout << deltamax << endl;
}

void LidDrivenCavity::solvePoisson()
{
    // TODO: Collect and deserialize into the root process

    // Fill up those buffers!
    filltxbuf(v_l);
    MPI_Allgatherv(txrxBuffer, (Nyl * Nxl), MPI_DOUBLE, v_s, txrxcounts, displs, MPI_DOUBLE, grid);
    deserialize(v_s, v);

    // Get inner sections...
    double *b = new double[(Ny - 2) * (Nx - 2)];
    double *s_result = new double[(Ny - 2) * (Nx - 2)];

    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            b[(j - 1) + (i - 1) * (Ny - 2)] = v[j + i * Ny];
        }
    }

    pcgd(rank, P, dx2, dy2, (Nx - 2) * (Ny - 2), Ny - 2, s_result, b);

    // Put s back where it belongs

    for (int j = 1; j < Ny - 1; j++)
    {
        for (int i = 1; i < Nx - 1; i++)
        {
            s[j + i * Ny] = s_result[(j - 1) + (i - 1) * (Ny - 2)];
        }
    }

    // Now re-scatter s by first serialising it
    serialize(s_s, s);
    MPI_Scatterv(s_s, txrxcounts, displs, MPI_DOUBLE, txrxBuffer, (Nyl * Nxl), MPI_DOUBLE, 0, grid);
    emptytxbuf(s_l);

    // Calculate central v and put it into a buffer
}

void LidDrivenCavity::filltxbuf(double *source)
{
    for (int j = 0; j < Nyl; j++)
    {
        for (int i = 0; i < Nxl; i++)
        {
            txrxBuffer[j + i * Nyl] = source[(j + 1) + (i + 1) * Nylh];
        }
    }
}

void LidDrivenCavity::emptytxbuf(double *dest)
{
    for (int j = 0; j < Nyl; j++)
    {
        for (int i = 0; i < Nxl; i++)
        {
            dest[(j + 1) + (i + 1) * Nylh] = txrxBuffer[j + i * Nyl];
        }
    }
}

void LidDrivenCavity::deserialize(double *v_s, double *v)
{
    // Organise this back into matrix v from v_s...
    int cursor = 0;

    for (int i = 0; i < Px; i++)
    {
        for (int j = 0; j < Py; j++)
        {
            // For each cell.. Figure out its coordinates

            int cell_y = (i == Px - 1) ? (Nxb + rem_x) : Nxb;
            int cell_x = (j == Py - 1) ? (Nyb + rem_y) : Nyb;

            int cursor_serial = displs[j + i * Py];

            // x-y position of origin of submatrix
            int cursor_matx = j * (Nxb);
            int cursor_maty = i * (Nyb);

            // cout << cursor_matx << " " << cursor_maty << endl;

            // k <- local i
            // l <- local j
            for (int k = 0; k < cell_x; k++)
            {
                for (int l = 0; l < cell_y; l++)
                {

                    v[(cursor_maty + l) + (cursor_matx + k) * Ny] = v_s[cursor_serial + (l + k * cell_y)];
                }
            }
        }
    }
}

void LidDrivenCavity::serialize(double *v_s, double *v)
{
    // Organise this back into matrix v from v_s...
    int cursor = 0;

    for (int i = 0; i < Px; i++)
    {
        for (int j = 0; j < Py; j++)
        {
            // For each cell.. Figure out its coordinates

            int cell_y = (i == Px - 1) ? (Nxb + rem_x) : Nxb;
            int cell_x = (j == Py - 1) ? (Nyb + rem_y) : Nyb;

            int cursor_serial = displs[j + i * Py];

            // x-y position of origin of submatrix
            int cursor_matx = j * (Nxb);
            int cursor_maty = i * (Nyb);

            // cout << cursor_matx << " " << cursor_maty << endl;

            // k <- local i
            // l <- local j
            for (int k = 0; k < cell_x; k++)
            {
                for (int l = 0; l < cell_y; l++)
                {
                    v_s[cursor_serial + (l + k * cell_y)] = v[(cursor_maty + l) + (cursor_matx + k) * Ny];
                }
            }
        }
    }
}

// Send a message to to each of my neighbours,
// sending my boundaries and asking for my new boundaries
void LidDrivenCavity::updateHalo()
{
    updatevHalo();
    updatesHalo();
}

void LidDrivenCavity::updatevHalo()
{
    // Fill in my transmit buffers
    for (int i = 0; i < Nxl; i++)
    {
        txtopBuffer[i] = v_l[(Nyl - 1 + 1) + (i + 1) * Nylh];
        txbottomBuffer[i] = v_l[(0 + 1) + (i + 1) * Nylh];
    }

    for (int j = 0; j < Nyl; j++)
    {
        txleftBuffer[j] = v_l[(j + 1) + (0 + 1) * Nylh];
        txrightBuffer[j] = v_l[(j + 1) + (Nxl - 1 + 1) * Nylh];
    }

    txrx();

    // Put those received buffers into my own local matrices
    for (int i = 0; i < Nxl; i++)
    {
        v_l[(Nyl - 1 + 2) + (i + 1) * Nylh] = rxtopBuffer[i];
        v_l[(0) + (i + 1) * Nylh] = rxbottomBuffer[i];
    }

    for (int j = 0; j < Nyl; j++)
    {
        v_l[(j + 1) + (0) * Nylh] = rxleftBuffer[j];
        v_l[(j + 1) + (Nxl - 1 + 2) * Nylh] = txrightBuffer[j];
    }
}

void LidDrivenCavity::updatesHalo()
{
    // Fill in my transmit buffers
    for (int i = 0; i < Nxl; i++)
    {
        txtopBuffer[i] = s_l[(Nyl - 1 + 1) + (i + 1) * Nylh];
        txbottomBuffer[i] = s_l[(0 + 1) + (i + 1) * Nylh];
    }

    for (int j = 0; j < Nyl; j++)
    {
        txleftBuffer[j] = s_l[(j + 1) + (0 + 1) * Nylh];
        txrightBuffer[j] = s_l[(j + 1) + (Nxl - 1 + 1) * Nylh];
    }

    txrx();

    // Put those received buffers into my own local matrices
    for (int i = 0; i < Nxl; i++)
    {
        s_l[(Nyl - 1 + 2) + (i + 1) * Nylh] = rxtopBuffer[i];
        s_l[(0) + (i + 1) * Nylh] = rxbottomBuffer[i];
    }

    for (int j = 0; j < Nyl; j++)
    {
        s_l[(j + 1) + (0) * Nylh] = rxleftBuffer[j];
        s_l[(j + 1) + (Nxl - 1 + 2) * Nylh] = txrightBuffer[j];
    }
}

void LidDrivenCavity::txrx()
{
    // Unfortunately can't condense these as Recv is blocking
    // Tx/Rx
    if (neighbour_up != -1)
    {
        // Send non-blocking
        MPI_Send(txtopBuffer, Nxl, MPI_DOUBLE, neighbour_up, 0, grid);
    }

    // Send to lower neighbours
    if (neighbour_down != -1)
    {
        // Send non-blocking
        MPI_Send(txbottomBuffer, Nxl, MPI_DOUBLE, neighbour_down, 0, grid);
    }

    if (neighbour_left != -1)
    {
        // Send non-blocking
        MPI_Send(txleftBuffer, Nyl, MPI_DOUBLE, neighbour_left, 0, grid);
    }

    if (neighbour_right != -1)
    {
        // Send non-blocking
        MPI_Send(txrightBuffer, Nyl, MPI_DOUBLE, neighbour_right, 0, grid);
    }

    if (neighbour_up != -1)
    {
        MPI_Recv(rxtopBuffer, Nxl, MPI_DOUBLE, neighbour_up, 0, grid, MPI_STATUS_IGNORE);
    }

    if (neighbour_down != -1)
    {
        MPI_Recv(rxbottomBuffer, Nxl, MPI_DOUBLE, neighbour_down, 0, grid, MPI_STATUS_IGNORE);
    }

    if (neighbour_left != -1)
    {
        MPI_Recv(rxleftBuffer, Nyl, MPI_DOUBLE, neighbour_left, 0, grid, MPI_STATUS_IGNORE);
    }

    if (neighbour_right != -1)
    {
        MPI_Recv(rxrightBuffer, Nyl, MPI_DOUBLE, neighbour_right, 0, grid, MPI_STATUS_IGNORE);
    }
}