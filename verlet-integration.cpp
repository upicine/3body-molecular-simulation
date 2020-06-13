#include <cmath>
#include <mpi.h>
#include <cassert>

#include "verlet-integration.h"
#include "utils.h"
#include "embedded-algorithm.h"
#include "particle-buffer.h"


static double calcDistance(Particle &i, Particle &j) {
    double dist = std::sqrt(std::pow(i.x.coor - j.x.coor, 2)
                            + std::pow(i.y.coor - j.y.coor, 2)
                            + std::pow(i.z.coor - j.z.coor, 2));
    return std::max(dist, MIN_ABS);
}

static double calcPotential1D(Particle1D &i1D, Particle &i, Particle &j, Particle &k) {
    double rij = calcDistance(i, j);
    double rik = calcDistance(i, k);
    double rkj = calcDistance(k, j);

    return (1 / std::pow(rij * rik * rkj, 3))
            + ((3
                * (-std::pow(rij, 2) + std::pow(rik, 2) + std::pow(rkj, 2))
                * (std::pow(rij, 2) - std::pow(rik, 2) + std::pow(rkj, 2))
                * (std::pow(rij, 2) + std::pow(rik, 2) - std::pow(rkj, 2)))
               / (8 * std::pow(rij * rik * rkj, 5)));
}

static void calcDerivative1D(Particle1D &i1D, Particle &i, Particle &j, Particle &k, bool twice) {
    double x = i1D.coor;
    double h = std::abs(x) < MIN_ABS ?
               EPS * MIN_ABS :
               EPS * x;
    i1D.coor = x + h;
    double fxph = calcPotential1D(i1D, i, j, k);
    i1D.coor = x - h;
    double fxmh = calcPotential1D(i1D, i, j, k);
    i1D.coor = x;

    i1D.f += twice ? 2.0 * (fxph - fxmh) : (fxph - fxmh);
}

static void calcDerivative(Particle &i, Particle &j, Particle &k, bool twice=false) {
    calcDerivative1D(i.x, i, j, k, twice);
    calcDerivative1D(i.y, i, j, k, twice);
    calcDerivative1D(i.z, i, j, k, twice);
}

void computeForce(ParticleBuff &b1, ParticleBuff &b2, ParticleBuff &b3) {
    assert(b1.buf[0].id==startIndex(b1.owner, b1.p, b1.n));
    assert(b2.buf[0].id==startIndex(b2.owner, b2.p, b2.n));
    assert(b3.buf[0].id==startIndex(b3.owner, b3.p, b3.n));

    for (int i = 0; i < b1.buf_sz; i++) {
        for (int j = 0; j < b2.buf_sz; j++) {
            for (int k = 0; k < b3.buf_sz; k++) {
                if (b1.owner == b2.owner && b2.owner == b3.owner) {
                    if (i != j && j < k && i != k)
                        calcDerivative(b1.buf[i], b2.buf[j], b3.buf[k], true);
                } else if (b1.owner == b2.owner) {
                    if (i != j) {
                        calcDerivative(b1.buf[i], b2.buf[j], b3.buf[k], true);
                    }
                    if (i < j) {
                        calcDerivative(b3.buf[k], b1.buf[i], b2.buf[j], true);
                    }
                } else if (b2.owner == b3.owner) {
                    if (j != k) {
                        calcDerivative(b2.buf[j], b1.buf[i], b3.buf[k], true);
                    }
                    if (j < k) {
                        calcDerivative(b1.buf[i], b2.buf[j], b3.buf[k], true);
                    }
                } else if (b1.owner == b3.owner) {
                    if (i != k) {
                        calcDerivative(b1.buf[i], b2.buf[j], b3.buf[k], true);
                    }
                    if (i < k) {
                        calcDerivative(b2.buf[j], b1.buf[i], b3.buf[k], true);
                    }
                } else {
                    calcDerivative(b2.buf[j], b1.buf[i], b3.buf[k], true);
                    calcDerivative(b1.buf[i], b2.buf[j], b3.buf[k], true);
                    calcDerivative(b3.buf[k], b2.buf[j], b1.buf[i], true);
                }
            }
        }
    }
}


static void sumForce1D(Particle1D &i, Particle1D &j, Particle1D &k) {
    double x = i.coor;
    double h = std::abs(x) < MIN_ABS ?
               EPS * MIN_ABS :
               EPS * x;
    volatile double dx = (x + h) - (x - h);
    i.f += j.f + k.f;
    i.f /= dx;
}


void sumForces(ParticleBuff* b) {
    for (int i = 0; i < b->buf_sz; i++) {
        sumForce1D(b[0].buf[i].x, b[1].buf[i].x, b[2].buf[i].x);
        sumForce1D(b[0].buf[i].y, b[1].buf[i].y, b[2].buf[i].y);
        sumForce1D(b[0].buf[i].z, b[1].buf[i].z, b[2].buf[i].z);
    }
}


void calcStartingAcc1D(Particle1D &p) {
    p.a = -(p.f / M);
    p.f = 0;
}


void calcStartingAcc(Particle *particles, int particles_sz) {
    for (int i = 0; i < particles_sz; i++) {
        calcStartingAcc1D(particles[i].x);
        calcStartingAcc1D(particles[i].y);
        calcStartingAcc1D(particles[i].z);
    }
}


void calcNewCoor1D(Particle1D &p, double dt) {
    p.coor = p.coor + p.v * dt + p.a * dt * dt / 2;
}


void calcNewCoor(Particle *particles, int particles_sz, double dt) {
    for (int i = 0; i < particles_sz; i++) {
        calcNewCoor1D(particles[i].x, dt);
        calcNewCoor1D(particles[i].y, dt);
        calcNewCoor1D(particles[i].z, dt);
    }
}


void calcNewAccAndV1D(Particle1D &p, double dt) {
    double da = -p.f / M;
    p.v = p.v + (p.a + da) * dt / 2;
    p.a = da;
    p.f = 0;
}


void calcNewAccAndV(Particle *particles, int particles_sz, double dt) {
    for (int i = 0; i < particles_sz; i++) {
        calcNewAccAndV1D(particles[i].x, dt);
        calcNewAccAndV1D(particles[i].y, dt);
        calcNewAccAndV1D(particles[i].z, dt);
    }
}