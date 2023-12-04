/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-07 21:41:18
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: test file for ParticleStruct
 */

#ifndef PARTICLE_STRUCT_TEST_CPP_
#define PARTICLE_STRUCT_TEST_CPP_

#include <iostream>
#include "Particle/Particle.hpp"

int main() {
    Particle::Particle p1;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Particle p1;" << std::endl;
    std::cout << p1 << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Particle p2(3);" << std::endl;
    Particle::Particle p2(3);
    std::cout << p2 << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Particle p3(3, Particle::ParticleType::Wall);" << std::endl;
    Particle::Particle p3(3, Particle::ParticleType::Wall);
    std::cout << p3 << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Particle p4 = p3;" << std::endl;
    Particle::Particle p4;
    p4 = p3;
    std::cout << p4 << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "ParticleGroup pg(3, 3, Particle::ParticleType::Water);" << std::endl;
    Particle::ParticleGroup pg(3, 3, Particle::ParticleType::Water);
    for (int i = 0; i < pg.n_particles_; ++i) {
        pg[i].a_vec_ = i;
    }
    std::cout << pg << std::endl;
    return 0;
}

#endif // PARTICLE_STRUCT_TEST_CPP_