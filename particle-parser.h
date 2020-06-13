#ifndef MPI_HPC_PARTICLE_PARSER_H
#define MPI_HPC_PARTICLE_PARSER_H

#include "utils.h"

int parseParticles(std::string &filename, Particle **particles);

Particle *scatterParticles(Particle **particles, int n, int p, int rank,
                           int *send_counts, int *displacement);

int *generateBuffSizes(int p, int n);

int *generateDisplacement(int *send_counts, int p);

void gatherParticles(Particle *particles, Particle *my_particles,
                     int rank, int *recv_counts, int *displacement);

#endif //MPI_HPC_PARTICLE_PARSER_H
