//
// Created by michal on 11.06.20.
//

#ifndef MPI_HPC_VERLET_INTEGRATION_H
#define MPI_HPC_VERLET_INTEGRATION_H

#include "embedded-algorithm.h"

void computeForce(ParticleBuff &b1, ParticleBuff &b2, ParticleBuff &b3);

void sumForces(ParticleBuff* b);

#endif //MPI_HPC_VERLET_INTEGRATION_H
