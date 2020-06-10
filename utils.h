#ifndef MPI_HPC_UTILS_H
#define MPI_HPC_UTILS_H

#include <algorithm>

struct __attribute__((packed)) Particle1D {
    double coor;
    double v;
    double f = 0;
    double a = 0;
};


struct __attribute__((packed)) Particle {
    Particle1D x;
    Particle1D y;
    Particle1D z;
};


int inline startIndex(int rank, int p, int n) {
    return rank * (n / p) + std::min(rank, n % p);
}


int inline nextRank(int rank, int p) {
    return rank == 0 ? p - 1 : rank - 1;
}


int inline prevRank(int rank, int p) {
    return (rank + 1) % p;
}


void printParticles(Particle *particles, size_t particles_sz);


#endif //MPI_HPC_UTILS_H
