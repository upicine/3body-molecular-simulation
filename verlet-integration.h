#ifndef MPI_HPC_VERLET_INTEGRATION_H
#define MPI_HPC_VERLET_INTEGRATION_H

#include "embedded-algorithm.h"
#include "particle-buffer.h"

void computeForce(ParticleBuff &b1, ParticleBuff &b2, ParticleBuff &b3);

void sumForces(ParticleBuff* b);

void calcStartingAcc(Particle *particles, int particles_sz);

void calcNewCoor(Particle *particles, int particles_sz, double dt);

void calcNewAccAndV(Particle *particles, int particles_sz, double dt);

void computeForceSeq(Particle *particles, int particles_sz);

#endif //MPI_HPC_VERLET_INTEGRATION_H
