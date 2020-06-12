#include <vector>
#include <fstream>
#include <sstream>
#include <numeric>
#include <mpi.h>

#include "particle-parser.h"

static int parseParticles(const char* filename, Particle **particles) {
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

    return static_cast<int>(particles_vec.size());
}

void generateBuffSizes(int *send_counts, int p, int n) {
    for (int i = 0; i < p; i++) {
        send_counts[i] = bufferSize(i, p, n) * PARTICLE_SIZE;
    }
}

void generateDisplacement(int *displacement, int *send_counts, int p) {
    displacement[0] = 0;
    std::partial_sum(send_counts, send_counts + p - 1, displacement + 1);
}

void scatterParticles(Particle **particles, Particle **my_particles, int n,
                      int p, int rank) {
    int *send_counts = new int[p];
    int *displacement = new int[p];

    generateBuffSizes(send_counts, p, n);
    generateDisplacement(displacement, send_counts, p);

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

    *my_particles = (Particle*)(recv_buff);
}

void gatherParticles
