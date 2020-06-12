#ifndef MPI_HPC_PARTICLE_PARSER_H
#define MPI_HPC_PARTICLE_PARSER_H

#include "utils.h"

int parseParticles(const char* filename, Particle **particles);

void scatterParticles(Particle **particles, Particle **my_particles, int n,
                      int p, int rank);

void generateBuffSizes(int *send_counts, int p, int n);

void generateDisplacement(int *displacement, int *send_counts, int p);

#endif //MPI_HPC_PARTICLE_PARSER_H
