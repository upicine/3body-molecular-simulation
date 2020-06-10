#include <vector>
#include <fstream>
#include <sstream>

#include "particle-parser.h"

size_t parseParticles(const char* filename, Particle **particles) {
    std::vector<Particle> particles_vec;
    std::ifstream particles_file(filename);

    for (std::string line; std::getline(particles_file, line); ) {
        Particle p;
        std::stringstream ss(line);

        ss >> p.x.coor >> p.y.coor >> p.z.coor;
        ss >> p.x.v >> p.y.v >> p.z.v;

        particles_vec.push_back(p);
    }

    *particles = new Particle[particles_vec.size()];
    std::copy(particles_vec.begin(), particles_vec.end(), *particles);

    return particles_vec.size();
}

