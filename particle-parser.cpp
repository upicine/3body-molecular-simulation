#include <vector>
#include <fstream>
#include <sstream>
#include <numeric>
#include <mpi.h>

#include "particle-parser.h"

int parseParticles(std::string &filename, Particle **particles) {
    std::vector<Particle> particles_vec;
    std::ifstream particles_file(filename);
    double id_counter = 0.0;

    for (std::string line; std::getline(particles_file, line); ) {
        Particle p;
        std::stringstream ss(line);

        ss >> p.x.coor >> p.y.coor >> p.z.coor;
        ss >> p.x.v >> p.y.v >> p.z.v;
        p.id = id_counter;
        id_counter += 1.0;

        particles_vec.push_back(p);
    }

    *particles = new Particle[particles_vec.size()];
    std::copy(particles_vec.begin(), particles_vec.end(), *particles);
    particles_file.close();

    return static_cast<int>(particles_vec.size());
}

int* generateBuffSizes(int p, int n) {
    int *send_counts = new int[p];

    for (int i = 0; i < p; i++) {
        send_counts[i] = bufferSize(i, p, n) * PARTICLE_SIZE;
    }

    return send_counts;
}

int* generateDisplacement(int *send_counts, int p) {
    int *displacement = new int[p];

    displacement[0] = 0;
    std::partial_sum(send_counts, send_counts + p - 1, displacement + 1);

    return displacement;
}

Particle* scatterParticles(Particle **particles, int n, int p, int rank,
                           int* send_counts, int* displacement) {
    int recv_buff_sz = bufferSize(rank, p, n) * PARTICLE_SIZE;
    double *recv_buff = new double[recv_buff_sz];

    double *part_as_double = (double *) (*particles);
    MPI_Scatterv(part_as_double,
                 send_counts,
                 displacement,
                 MPI_DOUBLE,
                 recv_buff,
                 recv_buff_sz,
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);

    return (Particle*)(recv_buff);
}

void gatherParticles(Particle *particles, Particle *my_particles,
                     int rank, int* recv_counts, int* displacement) {
    double *recv_buf = (double*)particles;
    double *send_buf = (double*)my_particles;
    int send_sz = recv_counts[rank];

    MPI_Gatherv(send_buf,
                send_sz,
                MPI_DOUBLE,
                recv_buf,
                recv_counts,
                displacement,
                MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

}
