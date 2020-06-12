#include "particle-buffer.h"

ParticleBuff::ParticleBuff(int owner, int p, int n)
        : owner(owner), p(p), n(n) {
    buf = new Particle[bufferSize(0, p, n)];
    d_recv_buf = new double[bufferSize(0, p, n) * PARTICLE_SIZE];
    buf_sz = bufferSize(owner, p, n);
    d_buf = (double*)buf;
    d_buf_sz = buf_sz * PARTICLE_SIZE;
}

ParticleBuff::ParticleBuff(Particle *my_buf, int owner, int p, int n)
        : owner(owner), p(p), n(n) {
    buf = new Particle[bufferSize(0, p, n)];
    d_recv_buf = new double[bufferSize(0, p, n) * PARTICLE_SIZE];
    buf_sz = bufferSize(owner, p, n);
    d_buf = (double*)buf;
    d_buf_sz = buf_sz * PARTICLE_SIZE;
    std::copy(my_buf, my_buf + buf_sz, buf);
}

ParticleBuff::~ParticleBuff() {
    delete[] buf;
    delete[] d_recv_buf;
}

void ParticleBuff::setNext() {
    owner = nextRank(owner, p);
    buf_sz = bufferSize(owner, p, n);
    d_buf_sz = buf_sz * PARTICLE_SIZE;
}

void ParticleBuff::setPrev() {
    owner = prevRank(owner, p);
    buf_sz = bufferSize(owner, p, n);
    d_buf_sz = buf_sz * PARTICLE_SIZE;
}

void ParticleBuff::setOwner(int owner_) {
    owner = owner_;
    buf_sz = bufferSize(owner, p, n);
    d_buf_sz = buf_sz * PARTICLE_SIZE;
}

void ParticleBuff::switchRecvBuf() {
    double *tmp = d_buf;
    d_buf = d_recv_buf;
    buf = (Particle*)d_recv_buf;
    d_recv_buf = tmp;
}

