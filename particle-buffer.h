//
// Created by michal on 12.06.20.
//

#ifndef MPI_HPC_PARTICLE_BUFFER_H
#define MPI_HPC_PARTICLE_BUFFER_H

#include "utils.h"

class ParticleBuff {
public:
    Particle *buf;
    double *d_buf;
    double *d_recv_buf;
    int d_buf_sz;
    int buf_sz;
    int owner;
    int p;
    int n;

    ParticleBuff(int owner, int p, int n);

    ParticleBuff(Particle *my_buf, int owner, int p, int n);

    ~ParticleBuff();

    void setNext();

    void setPrev();

    void setOwner(int owner_);

    void switchRecvBuf();

    int inline getNext() {
        return nextRank(owner, p);
    }

    int inline getPrev() {
        return prevRank(owner, p);
    }

};

#endif //MPI_HPC_PARTICLE_BUFFER_H
