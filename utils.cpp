#include "utils.h"

#include <iostream>
#include <iterator>

static void printParticle(Particle p) {
    std::cout << p.x.coor << " " << p.y.coor << " " << p.z.coor << " "
              << p.x.v << " " << p.y.v << " " << p.z.v << std::endl;
}

void printParticles(Particle *particles, int particles_sz) {
    for (int i = 0; i < particles_sz; i++) {
        printParticle(particles[i]);
    }
}