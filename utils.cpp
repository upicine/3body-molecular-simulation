#include <iostream>
#include <iterator>
#include <cstring>
#include <fstream>
#include <iomanip>

#include "utils.h"

static void printParticle(Particle p) {
    std::cout << p.x.coor << " " << p.y.coor << " " << p.z.coor << " "
              << p.x.v << " " << p.y.v << " " << p.z.v << std::endl;
}

void printParticles(Particle *particles, int particles_sz) {
    for (int i = 0; i < particles_sz; i++) {
        printParticle(particles[i]);
    }
}

void parseArgs(int argc, char** argv, double &dt, int &steps,
               std::string &input, std::string &output, bool &verbose) {
    input = argv[1];
    output = argv[2];
    steps = std::atoi(argv[3]);
    dt = std::atof(argv[4]);
    verbose = argc == 6 ? std::strcmp(argv[5], "-v") == 0 : false;
}

void saveResults(Particle *particles, int particle_sz, std::string &filename,
                 int step) {
    std::string output_filename = filename + "_" + std::to_string(step) + ".txt";
    std::ofstream output_file(output_filename);
    output_file << std::setprecision(16);

    for (int i = 0; i < particle_sz; i++) {
        output_file << particles[i].x.coor << " "
                    << particles[i].y.coor << " "
                    << particles[i].z.coor << " "
                    << particles[i].x.v << " "
                    << particles[i].y.v << " "
                    << particles[i].z.v << std::endl;
    }
}