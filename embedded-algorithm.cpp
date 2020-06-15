#include <mpi.h>

#include "embedded-algorithm.h"
#include "verlet-integration.h"
#include "particle-buffer.h"

const int START_BUFF_TAG = 0;
const int SHIFT_BUFF_TAG = 1;
const int SHIFT_BUFF_TAG_2 = 2;
const int RECV_BUFF_TAG = 3;


void shiftRight(ParticleBuff &pb, int tag, int rank, int p) {
    MPI_Request request[2];
    MPI_Status status[2];
    int next = nextRank(rank, p);
    int prev = prevRank(rank, p);

    MPI_Isend(pb.d_buf, pb.d_buf_sz, MPI_DOUBLE, next, tag, MPI_COMM_WORLD, &request[0]);
    pb.setPrev();
    MPI_Irecv(pb.d_recv_buf, pb.d_buf_sz, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, &request[1]);
    MPI_Waitall(2, request, status);

    pb.switchRecvBuf();
}

void sendAndRecvResults(ParticleBuff *b, int tag, int rank) {
    MPI_Request request[6];
    MPI_Status status[6];

    for (int i = 0, j = 0; i < 3; i++, j+= 2) {
        MPI_Isend(b[i].d_buf, b[i].d_buf_sz, MPI_DOUBLE, b[i].owner, tag, MPI_COMM_WORLD, &request[j]);
        b[i].setOwner(rank);
        MPI_Irecv(b[i].d_recv_buf, b[i].d_buf_sz, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &request[j + 1]);
    }
    MPI_Waitall(6, request, status);

    for (int i = 0; i < 3; i++) {
        b[i].switchRecvBuf();
    }
}

void embeddedAlgorithm(Particle *particles, int rank, int p, int n) {
    ParticleBuff b[] = {
            {prevRank(rank, p), p, n},
            {particles, rank, p, n},
            {nextRank(rank, p), p, n}
    };

    int i = 0;

    MPI_Request request[4];
    MPI_Status status[4];

    MPI_Isend(b[1].d_buf, b[1].d_buf_sz, MPI_DOUBLE, b[1].getPrev(), START_BUFF_TAG, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(b[1].d_buf, b[1].d_buf_sz, MPI_DOUBLE, b[1].getNext(), START_BUFF_TAG, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(b[2].d_buf, b[2].d_buf_sz, MPI_DOUBLE, b[2].owner, START_BUFF_TAG, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(b[0].d_buf, b[0].d_buf_sz, MPI_DOUBLE, b[0].owner, START_BUFF_TAG, MPI_COMM_WORLD, &request[3]);
    MPI_Waitall(4, request, status);

    for (int s = p - 3; s >= 0; s -= 3) {
        for (int j = 0; j < s; j++) {
            if (j != 0 || s != p - 3) {
                shiftRight(b[i], SHIFT_BUFF_TAG, rank, p);
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
        i = i == 0 ? 2 : i - 1;

        shiftRight(b[i], SHIFT_BUFF_TAG_2, rank, p);

        if ((rank / (p / 3)) == 0) {
            computeForce(b[0], b[1], b[2]);
        }

    }

    sendAndRecvResults(b, RECV_BUFF_TAG, rank);
    sumForces(b);
    std::copy(b[0].buf, b[0].buf + b[0].buf_sz, particles);
}