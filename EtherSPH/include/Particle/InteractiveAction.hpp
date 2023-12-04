/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-08 15:59:23
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: header file for InteractiveAction
 */

#ifndef INTERACTIVE_ACTION_HPP_
#define INTERACTIVE_ACTION_HPP_

#include "ConstVariable/ConstPhysics.hpp"
#include "Container/Container.hpp"
#include "Particle/ParticleStruct.hpp"

using namespace KernelFunction::Wendland;

namespace Particle {

void balanceMass(Particle& p_i, Particle& p_j, const double h, const int dim);
void balanceMass(ParticleGroup& particle_group, const double h, const int dim);

void internalForce(Particle& p_i, Particle& p_j, const double h, const int dim);
void internalForce(ParticleGroup& particle_group, const double h, const int dim);

void wallForce(Particle& p_i, const Container::DoubleArray& wall_vec, const double h, const int dim);
void wallForce(ParticleGroup& particle_group, ParticleGroup& wall_particle_group, const double h, const int dim);

};

#endif // INTERACTIVE_ACTION_HPP_