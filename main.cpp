#include <iostream>
#include <mpi.h>
#include <iomanip>

#include "utils.h"
#include "particle-parser.h"
#include "embedded-algorithm.h"
#include "verlet-integration.h"

int main(int argc, char *argv[]) {
    int num_processes, rank, steps;
    double dt, start_time;
    double par_exec = 0, seq_exec = 0;
    bool verbose = false;
    std::string input, output;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    if (rank == 0 && BENCHMARK) {
        start_time = MPI_Wtime();
    }

    parseArgs(argc, argv, dt, steps, input, output, verbose);

    std::cout << std::setprecision(16);

    Particle *particles, *my_particles;
    int particles_sz, my_particles_sz;

    if (rank == 0) {
        particles_sz = parseParticles(input, &particles);
    }

    if (num_processes == 1) {
        if (BENCHMARK) {
            seq_exec += MPI_Wtime() - start_time;
            start_time = MPI_Wtime();
        }

        for (int i = 0; i < steps; i++) {
            if (i == 0) {
                computeForceSeq(particles, particles_sz);
                calcStartingAcc(particles, particles_sz);
            }

            calcNewCoor(particles, particles_sz, dt);
            computeForceSeq(particles, particles_sz);
            calcNewAccAndV(particles, particles_sz, dt);

            if (verbose) {
                saveResults(particles, particles_sz, output, i + 1);
            }
        }
    } else {
        if (rank == 0 && BENCHMARK) {
            seq_exec += MPI_Wtime() - start_time;
            start_time = MPI_Wtime();
        }


        MPI_Bcast(&particles_sz, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int *send_counts = generateBuffSizes(num_processes, particles_sz);
        int *displacement = generateDisplacement(send_counts, particles_sz);

        my_particles = scatterParticles(&particles, particles_sz, num_processes,
                                        rank, send_counts, displacement);
        my_particles_sz = bufferSize(rank, num_processes, particles_sz);

        if (rank == 0 && BENCHMARK) {
            seq_exec += MPI_Wtime() - start_time;
            start_time = MPI_Wtime();
        }

        for (int i = 0; i < steps; i++) {
            if (i == 0) {
                embeddedAlgorithm(my_particles, rank, num_processes, particles_sz);
                calcStartingAcc(my_particles, my_particles_sz);
            }

            calcNewCoor(my_particles, my_particles_sz, dt);
            embeddedAlgorithm(my_particles, rank, num_processes, particles_sz);
            calcNewAccAndV(my_particles, my_particles_sz, dt);

            if (verbose) {
                gatherParticles(particles, my_particles, rank,
                                send_counts, displacement);
                if (rank == 0) {
                    saveResults(particles, particles_sz, output, i + 1);
                }
            }
        }

        gatherParticles(particles, my_particles, rank, send_counts, displacement);
    }

    if (rank == 0 && BENCHMARK) {
        par_exec += MPI_Wtime() - start_time;
        start_time = MPI_Wtime();
    }

    if (rank == 0) {
        saveResults(particles, particles_sz, output, steps);
    }

    if (rank == 0 && BENCHMARK) {
        seq_exec += MPI_Wtime() - start_time;
        std::cout << "seq_exec=" << seq_exec << " " << "par_exec=" << par_exec << " seq+par=" << seq_exec + par_exec << std::endl;
    }

    MPI_Finalize();

    return 0;
}


