/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-07 21:26:25
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: source file for SelfAction
 */

#ifndef SELF_ACTION_CPP_
#define SELF_ACTION_CPP_

#include "Particle/SelfAction.hpp"

namespace Particle {

void updateVelocity(Particle& particle, const double dt) {
    particle.v_vec_ += particle.a_vec_ * dt;
    return;
};

void updateVelocity(ParticleGroup& particle_group, const double dt) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        updateVelocity(particle_group[i], dt);
    }
    return;
};

void updateVelocity(Particle& particle, const double dt, const Container::DoubleArray& body_force_vec) {
    particle.a_vec_ += body_force_vec;
    particle.v_vec_ += particle.a_vec_ * dt;
    return;
};

void updateVelocity(ParticleGroup& particle_group, const double dt, const Container::DoubleArray& body_force_vec) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        updateVelocity(particle_group[i], dt, body_force_vec);
    }
    return;
};

void updatePosition(Particle& particle, const double dt) {
    particle.x_vec_ += (particle.v_vec_ - particle.a_vec_ * dt / 2) * dt;
    particle.a_vec_ = 0.;
    // usually, update position is the last step of a time step loop
    // thus we need to reset acceleration to 0
    return;
};

void updatePosition(ParticleGroup& particle_group, const double dt) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        updatePosition(particle_group[i], dt);
    }
    return;
};

void updateDensity(Particle& particle, const double dt) {
    particle.rho_ += particle.drho_ * dt;
    particle.drho_ = 0.;
    return;
};

void updateDensity(ParticleGroup& particle_group, const double dt) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        updateDensity(particle_group[i], dt);
    }
    return;
};

void updatePressure(Particle& particle, const double dt) {
    particle.rho_ += particle.drho_ * dt;
    particle.drho_ = 0.;
    particle.p_ = ConstVariable::Physics::Water::c_02 * 
                    (particle.rho_ - ConstVariable::Physics::Water::rho_0) +
                        ConstVariable::Physics::Water::p_0;
    return;
};

void updatePressure(ParticleGroup& particle_group, const double dt) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        updatePressure(particle_group[i], dt);
    }
    return;
};

};

#endif // SELF_ACTION_CPP_