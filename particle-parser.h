#ifndef MPI_HPC_PARTICLE_PARSER_H
#define MPI_HPC_PARTICLE_PARSER_H

#include "utils.h"

int parseParticles(const char* filename, Particle **particles);

void scatterParticles(Particle **particles, Particle **my_particles, int n,
                      int p, int rank);

#endif //MPI_HPC_PARTICLE_PARSER_H
