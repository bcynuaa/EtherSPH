/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-07 21:24:38
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: header file for SelfAction
 */

#ifndef SELF_ACTION_HPP_
#define SELF_ACTION_HPP_

#include "Particle/ParticleStruct.hpp"

namespace Particle {

void updateVelocity(Particle& particle, const double dt);
void updateVelocity(ParticleGroup& particle_group, const double dt);
void updateVelocity(Particle& particle, const double dt, const Container::DoubleArray& body_force_vec);
void updateVelocity(ParticleGroup& particle_group, const double dt, const Container::DoubleArray& body_force_vec);

void updatePosition(Particle& particle, const double dt);
void updatePosition(ParticleGroup& particle_group, const double dt);

void updateDensity(Particle& particle, const double dt);
void updateDensity(ParticleGroup& particle_group, const double dt);

void updatePressure(Particle& particle, const double dt);
void updatePressure(ParticleGroup& particle_group, const double dt);

};

#endif // SELF_ACTION_HPP_