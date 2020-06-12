//
// Created by michal on 11.06.20.
//

#ifndef MPI_HPC_EMBEDDED_ALGORITHM_H
#define MPI_HPC_EMBEDDED_ALGORITHM_H

#include "utils.h"

class ParticleBuff {
public:
    Particle *buf;
    double *d_buf;
    int d_buf_sz;
    int buf_sz;
    int owner;
    int p;
    int n;

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
        d_buf_sz = buf_sz * PARTICLE_SIZE;
    }

    void inline setPrev() {
        owner = prevRank(owner, p);
        buf_sz = bufferSize(owner, p, n);
        d_buf_sz = buf_sz * PARTICLE_SIZE;
    }

    void inline setOwner(int owner_) {
        owner = owner_;
        buf_sz = bufferSize(owner, p, n);
        d_buf_sz = buf_sz * PARTICLE_SIZE;
    }

    int inline getNext() {
        return nextRank(owner, p);
    }

    int inline getPrev() {
        return prevRank(owner, p);
    }

};

void embeddedAlgorithm(Particle *particles, int rank, int p, int n);

#endif //MPI_HPC_EMBEDDED_ALGORITHM_H
