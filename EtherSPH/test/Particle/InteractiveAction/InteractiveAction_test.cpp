/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-08 21:01:08
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: test file for InteractiveAction
 */

#ifndef INTERACTIVE_ACTION_TEST_CPP_
#define INTERACTIVE_ACTION_TEST_CPP_

#include <iostream>
#include "Container/Container.hpp"
#include "Particle/Particle.hpp"

int main(int argc, char const *argv[]) {
    Container::DoubleArray gravity_vec = {0., -9.8};
    const double dt = 0.01;
    const double h = 0.1;
    const double dy = 0.03;
    const int dim = 2;
    const int n_particles = 3;
    Particle::ParticleGroup pg(dim, n_particles, Particle::ParticleType::Water);
    for (int i = 0; i < n_particles; ++i) {
        pg[i].x_vec_ = {0., -i * dy};
        pg[i].v_vec_ = {i*1., 0.};
        pg[i].a_vec_ = {i*2., i*3.};
        pg[i].p_ = ConstVariable::Physics::Water::rho_0 * 
                    ConstVariable::Physics::gravity * i * dy + 
                    ConstVariable::Physics::Water::p_0;
        pg[i].rho_ = ConstVariable::Physics::Water::rho_0 +
                        (pg[i].p_ - ConstVariable::Physics::Water::p_0)
                            / ConstVariable::Physics::Water::c_02;
    }
    std::cout << pg << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "start a virtual loop" << std::endl;
    std::cout << "Particle::balanceMass(Particle& p_i, Particle& p_j, const double h, const int dim)" << std::endl;
    for (int i = 0; i < n_particles; i++) {
        for (int j = 0; j < i+1; j++) {
            Particle::balanceMass(pg[i], pg[j], h, dim);
        }
    }
    std::cout << "updateDensity(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updateDensity(pg, dt);
    std::cout << "updatePressure(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updatePressure(pg, dt);
    std::cout << "Particle::internalsForce(Particle& p_i, Particle& p_j, const double h, const int dim)" << std::endl;
    for (int i = 0; i < n_particles; i++) {
        for (int j = 0; j < i+1; j++) {
            Particle::internalForce(pg[i], pg[j], h, dim);
        }
    }
    std::cout << "updateVelocity(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updateVelocity(pg, dt, gravity_vec);
    std::cout << "updatePosition(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updatePosition(pg, dt);
    std::cout << "===================================================" << std::endl;
    std::cout << pg << std::endl;
    return 0;
}

#endif // INTERACTIVE_ACTION_TEST_CPP_