/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-07 21:58:08
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: test file for SelfAction
 */

#ifndef SELF_ACTION_TEST_CPP_
#define SELF_ACTION_TEST_CPP_

#include <iostream>
#include "ConstVariable/ConstVariable.hpp"
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
    std::cout << "updateDensity(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updateDensity(pg, dt);
    std::cout << pg << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "updatePressure(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updatePressure(pg, dt);
    std::cout << pg << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "updateVelocity(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updateVelocity(pg, dt, gravity_vec);
    std::cout << pg << std::endl;
    std::cout << "===================================================" << std::endl;
    std::cout << "updatePosition(ParticleGroup& particle_group, const double dt)" << std::endl;
    Particle::updatePosition(pg, dt);
    std::cout << pg << std::endl;
    std::cout << "For function for ParticleGroup relies on function for Particle" << std::endl;
    std::cout << "Tests only include function for ParticleGroup" << std::endl;
    return 0;
}

#endif // SELF_ACTION_TEST_CPP_