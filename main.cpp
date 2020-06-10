#include <iostream>
#include <mpi.h>

#include "utils.h"
#include "particle-parser.h"

int prev(int rank, int p) {
    return rank == 0 ? p - 1 : rank - 1;
}

int next(int rank, int p) {
//    std::cout << "--NEXT-- " << "rank=" << rank << " p=" << p << " result=" << (rank + 1) % p << std::endl;
    return (rank + 1) % p;
}

void compute(int rank, int a, int b, int c) {
    std::cout << "RANK=" << rank << " " << a << " " << b << " " << c << std::endl;
}

void calculateForces(int p, int rank) {
    int b[3];
    int i = 0;
    int left = prev(rank, p);
    int right = next(rank, p);
    b[1] = rank;

//    std::cout << "RANK=" << rank << " | " << "SEND TO " << left << std::endl;
    MPI_Send(&b[1], 1, MPI_INT, left, 0, MPI_COMM_WORLD);
//    std::cout << "RANK=" << rank << " | " << "SEND DONE " << std::endl;
//    std::cout << "RANK=" << rank << " | " << "SEND TO " << right << std::endl;
    MPI_Send(&b[1], 1, MPI_INT, right, 0, MPI_COMM_WORLD);
//    std::cout << "RANK=" << rank << " | " << "SEND DONE" << std::endl;
//    std::cout << "RANK=" << rank << " | " << "RECV FROM " << left << std::endl;
    MPI_Recv(&b[0], 1, MPI_INT, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    std::cout << "RANK=" << rank << " | " << "RECV DONE" << std::endl;
//    std::cout << "RANK=" << rank << " | " << "RECV FROM " << right << std::endl;
    MPI_Recv(&b[2], 1, MPI_INT, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    std::cout << "RANK=" << rank << " | " << "RECV DONE" << std::endl;

//    std::cout << "RANK=" << rank << " | " << b[0] << " " << b[1] << " " << b[2] << std::endl;

    for (int s = p - 3; s >= 0; s -= 3) {
        for (int j = 0; j < s; j++) {
            if (j != 0 || s != p - 3) {
                int i_rank = b[i];
//                std::cout << "LOOP RANK=" << rank << " i_rank=" << i_rank << std::endl;
//                std::cout << "LOOP RANK=" << rank << " | " << "SEND TO " << next(i_rank, p) << std::endl;
                MPI_Send(&b[i], 1, MPI_INT, next(i_rank, p), 1, MPI_COMM_WORLD);
//                std::cout << "LOOP RANK=" << rank << " | " << "SEND DONE " << std::endl;
//                std::cout << "LOOP RANK=" << rank << " | " << "RECV FROM " << prev(i_rank, p) << std::endl;
                MPI_Recv(&b[i], 1, MPI_INT, prev(i_rank, p), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                std::cout << "LOOP RANK=" << rank << " | " << "RECV DONE" << std::endl;
            } else {
                compute(rank, b[1], b[1], b[1]);
                compute(rank, b[1], b[1], b[2]);
                compute(rank, b[0], b[0], b[2]);
            }

            if (s == p - 3) {
                compute(rank, b[0], b[1], b[1]);
            }

            compute(rank, b[0], b[1], b[2]);
        }

        i = (i + 1) % 3;
    }

    if (p % 3 == 0) {
        i = (i - 1) % 3;
        int i_rank = b[i];

        MPI_Send(&b[i], 1, MPI_INT, next(i_rank, p), 2, MPI_COMM_WORLD);
        MPI_Recv(&b[i], 1, MPI_INT, prev(i_rank, p), 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if ((rank / (p / 3)) == 0) {
            compute(rank, b[0], b[1], b[2]);
        }

    }
}


int main(int argc, char *argv[]) {
    int num_processes, rank;

    const char* filename = "/home/michal/CLionProjects/mpi-hpc/test.in";

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    if (rank == 0) {
        Particle *particles;
        size_t particles_sz = parseParticles(filename, &particles);
        printParticles(particles, particles_sz);
    }

    MPI_Finalize();

    return 0;
}


