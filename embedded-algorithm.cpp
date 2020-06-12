#include <mpi.h>

#include "embedded-algorithm.h"
#include "verlet-integration.h"

static void compute(int rank, int a, int b, int c) {
    std::cout << "RANK " << rank << " | " << a << " " << b << " " << c << std::endl;
}

void shiftRight(ParticleBuff &pb, int tag) {
    MPI_Send(pb.d_buf, pb.d_buf_sz, MPI_DOUBLE, pb.getNext(), tag, MPI_COMM_WORLD);
    MPI_Recv(pb.d_buf, pb.d_buf_sz, MPI_DOUBLE, pb.getPrev(), tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    pb.setPrev();
}

void sendAndRecvResults(ParticleBuff *b, int tag, int rank) {
    MPI_Send(b[0].d_buf, b[0].d_buf_sz, MPI_DOUBLE, b[0].owner, tag, MPI_COMM_WORLD);
    MPI_Send(b[1].d_buf, b[1].d_buf_sz, MPI_DOUBLE, b[1].owner, tag, MPI_COMM_WORLD);
    MPI_Send(b[2].d_buf, b[2].d_buf_sz, MPI_DOUBLE, b[2].owner, tag, MPI_COMM_WORLD);
    b[0].setOwner(rank);
    b[1].setOwner(rank);
    b[2].setOwner(rank);
    MPI_Recv(b[0].d_buf, b[0].d_buf_sz, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(b[1].d_buf, b[1].d_buf_sz, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(b[2].d_buf, b[2].d_buf_sz, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}

void embeddedAlgorithm(Particle *particles, int rank, int p, int n) {
    ParticleBuff b[] = {
            {prevRank(rank, p), p, n},
            {particles, rank, p, n},
            {nextRank(rank, p), p, n}
    };

    int i = 0;

    MPI_Send(b[1].d_buf, b[1].d_buf_sz, MPI_DOUBLE, b[1].getPrev(), 0, MPI_COMM_WORLD);
    MPI_Send(b[1].d_buf, b[1].d_buf_sz, MPI_DOUBLE, b[1].getNext(), 0, MPI_COMM_WORLD);
    MPI_Recv(b[0].d_buf, b[0].d_buf_sz, MPI_DOUBLE, b[0].owner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(b[2].d_buf, b[2].d_buf_sz, MPI_DOUBLE, b[2].owner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int s = p - 3; s >= 0; s -= 3) {
        for (int j = 0; j < s; j++) {
            if (j != 0 || s != p - 3) {
                shiftRight(b[i], 1);
            } else {
                computeForce(b[1], b[1], b[1]);
                computeForce(b[1], b[1], b[2]);
                computeForce(b[0], b[0], b[2]);
            }

            if (s == p - 3) {
                computeForce(b[0], b[1], b[1]);
            }

            computeForce(b[0], b[1], b[2]);
        }

        i = (i + 1) % 3;
    }

    if (p % 3 == 0) {
        i = (i - 1) % 3;

        shiftRight(b[i], 2);

        if ((rank / (p / 3)) == 0) {
            computeForce(b[0], b[1], b[2]);
        }

    }

    sendAndRecvResults(b, 3, rank);
    sumForces(b);
    std::copy(b[0].buf, b[0].buf + b[0].buf_sz, particles);
}