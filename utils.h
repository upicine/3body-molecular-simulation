#ifndef MPI_HPC_UTILS_H
#define MPI_HPC_UTILS_H

#include <algorithm>
#include <iostream>
#include <iterator>

const int PARTICLE_SIZE = 13;
const double EPS = 4.69041575982343e-08;
const double MIN_ABS = 1e-10;
const double M = 1.0;

struct Particle1D {
    double coor;
    double v;
    double f = 0;
    double a = 0;
};

struct Particle {
    double id;
    Particle1D x;
    Particle1D y;
    Particle1D z;
};

int inline startIndex(int rank, int p, int n) {
    return rank * (n / p) + std::min(rank, n % p);
}

int inline bufferSize(int rank, int p, int n) {
    return n / p + (rank + 1 > (n % p) ? 0 : 1);
}

int inline prevRank(int rank, int p) {
    return rank == 0 ? p - 1 : rank - 1;
}

int inline nextRank(int rank, int p) {
    return (rank + 1) % p;
}

void printParticles(Particle *particles, int particles_sz);

template<typename T>
void printArray(T* arr, int size) {
    std::copy(arr, arr + size, std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

void parseArgs(int argc, char** argv, double &dt, int &steps,
               std::string &input, std::string &output, bool &verbose);

void saveResults(Particle *particles, int particle_sz, std::string &filename,
                 int step);


#endif //MPI_HPC_UTILS_H
