#include <iostream>
#include <mpi.h>
#include <iomanip>

#include "utils.h"
#include "particle-parser.h"
#include "embedded-algorithm.h"
#include "verlet-integration.h"

int main(int argc, char *argv[]) {
    int num_processes, rank;
    double dt = 0.5;

    const char* filename = "/home/michal/CLionProjects/mpi-hpc/test.in";

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    std::cout << std::setprecision(16);

    Particle *particles, *my_particles;
    int particles_sz, my_particles_sz;

    if (rank == 0) {
        particles_sz = parseParticles(filename, &particles);
    }

    MPI_Bcast(&particles_sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    scatterParticles(&particles, &my_particles, particles_sz, num_processes, rank);
    my_particles_sz = bufferSize(rank, num_processes, particles_sz);

    embeddedAlgorithm(my_particles, rank, num_processes, particles_sz);
    calcStartingAcc(my_particles, my_particles_sz);

    calcNewCoor(my_particles, my_particles_sz, dt);
    embeddedAlgorithm(my_particles, rank, num_processes, particles_sz);
    calcNewAccAndV(my_particles, my_particles_sz, dt);

    if (rank == 0) {
        printParticles(my_particles, my_particles_sz);
    }

    MPI_Finalize();

    return 0;
}


