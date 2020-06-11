#include <mpi.h>

#include "embedded-algorithm.h"

class ParticleBuff {
public:
    Particle *buf;
    double *d_buf;
    int d_buf_sz;
    int buf_sz;
    int owner;

    ParticleBuff(int owner, int p, int n) : owner(owner), p(p), n(n) {
        buf = new Particle[bufferSize(0, p, n)];
        buf_sz = bufferSize(owner, p, n);
        d_buf = (double*)buf;
        d_buf_sz = buf_sz * PARTICLE_SIZE;
    }

    ParticleBuff(Particle *my_buf, int owner, int p, int n) : owner(owner), p(p), n(n) {
        buf = new Particle[bufferSize(0, p, n)];
        buf_sz = bufferSize(owner, p, n);
        d_buf = (double*)buf;
        d_buf_sz = buf_sz * PARTICLE_SIZE;
        std::copy(my_buf, my_buf + buf_sz, buf);
    }

    ~ParticleBuff() {
        delete[] buf;
    }

    void inline setNext() {
        owner = nextRank(owner, p);
        buf_sz = bufferSize(owner, p, n);
    }

    void inline setPrev() {
        owner = prevRank(owner, p);
        buf_sz = bufferSize(owner, p, n);
    }

    int inline getNext() {
        return nextRank(owner, p);
    }

    int inline getPrev() {
        return prevRank(owner, p);
    }

private:
    int p;
    int n;

};

static void compute(int rank, int a, int b, int c) {
    std::cout << "RANK " << rank << " | " << a << " " << b << " " << c << std::endl;
}

void calculateForces(Particle *particles, int rank, int p, int n) {
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
                MPI_Send(b[i].d_buf, b[i].d_buf_sz, MPI_DOUBLE, b[i].getNext(), 1, MPI_COMM_WORLD);
                MPI_Recv(b[i].d_buf, b[i].d_buf_sz, MPI_DOUBLE, b[i].getPrev(), 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                b[i].setPrev();
            } else {
                compute(rank, b[1].owner, b[1].owner, b[1].owner);
                compute(rank, b[1].owner, b[1].owner, b[2].owner);
                compute(rank, b[0].owner, b[0].owner, b[2].owner);
            }

            if (s == p - 3) {
                compute(rank, b[0].owner, b[1].owner, b[1].owner);
            }

            compute(rank, b[0].owner, b[1].owner, b[2].owner);
        }

        i = (i + 1) % 3;
    }

    if (p % 3 == 0) {
        i = (i - 1) % 3;

        MPI_Send(b[i].d_buf, b[i].d_buf_sz, MPI_DOUBLE, b[i].getNext(), 2, MPI_COMM_WORLD);
        MPI_Recv(b[i].d_buf, b[i].d_buf_sz, MPI_DOUBLE, b[i].getPrev(), 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        b[i].setPrev();

        if ((rank / (p / 3)) == 0) {
            compute(rank, b[0].owner, b[1].owner, b[2].owner);
        }

    }
}