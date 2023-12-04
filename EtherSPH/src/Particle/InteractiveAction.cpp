/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-08 15:56:25
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: source file for InteractiveAction
 */

#ifndef INTERACTIVE_ACTION_CPP_
#define INTERACTIVE_ACTION_CPP_

#include "Particle/InteractiveAction.hpp"

using namespace KernelFunction::Wendland;

namespace Particle {

void balanceMass(Particle& p_i, Particle& p_j, const double h, const int dim) {
    const Container::DoubleArray delta_x_ij_vec = p_i.x_vec_ - p_j.x_vec_;
    const double r_ij = Container::norm(delta_x_ij_vec);
    if (r_ij > (radius_ratio * h)) {
        return;
    }
    else if (r_ij < ConstVariable::Math::epsilon) {
        // actually, the the 2 particles are the same particle
        return;
    }
    else {
        const double kernel_value = kernelValue(delta_x_ij_vec, h, dim);
        const Container::DoubleArray kernel_value_gradient = kernelValueGradient(delta_x_ij_vec, h, dim);
        const Container::DoubleArray delta_v_ij_vec = p_i.v_vec_ - p_j.v_vec_;
        const double drho = Container::dot(delta_v_ij_vec, delta_x_ij_vec);
        const double delta_rho_ij = p_i.rho_ - p_j.rho_;
        const double correction = 2.0 * ConstVariable::Physics::Water::epsilon_rho * 
                                    delta_rho_ij / r_ij * Container::norm(kernel_value_gradient);
        // according to the symmetrization of this part
        p_i.drho_ += p_j.mass_ * (drho + correction) * kernel_value;
        p_j.drho_ += p_i.mass_ * (drho - correction) * kernel_value;
        return;
    }
    return;
};

void balanceMass(ParticleGroup& particle_group, const double h, const int dim) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        for (int j = 0; j < i+1; ++j) {
            balanceMass(particle_group[i], particle_group[j], h, dim);
        }
    }
    return;
};

void internalForce(Particle& p_i, Particle& p_j, const double h, const int dim) {
    const Container::DoubleArray delta_x_ij_vec = p_i.x_vec_ - p_j.x_vec_;
    const double r_ij = Container::norm(delta_x_ij_vec);
    if (r_ij > (radius_ratio * h)) {
        return;
    }
    else if (r_ij < ConstVariable::Math::epsilon) {
        // actually, the the 2 particles are the same particle
        return;
    }
    else {
        const double dim_coefficient = 2.0 * (dim + 2.0);
        double coefficient = ConstVariable::Physics::Water::mu_0 * 
                                     dim_coefficient / p_i.rho_ / p_j.rho_;
        const Container::DoubleArray delta_v_ij_vec = p_i.v_vec_ - p_j.v_vec_;
        coefficient *= Container::dot(delta_v_ij_vec, delta_x_ij_vec) / 
                        (r_ij * r_ij + ConstVariable::Physics::avoid_singularity * h * h);
        const Container::DoubleArray kernel_value_gradient = kernelValueGradient(delta_x_ij_vec, h, dim);
        const Container::DoubleArray viscous_force_vec = coefficient * kernel_value_gradient;
        const Container::DoubleArray pressure_force_vec = (-p_i.p_/p_i.rho_/p_i.rho_  - 
                                                            p_j.p_/p_j.rho_/p_j.rho_) 
                                                                * kernel_value_gradient;
        p_i.a_vec_ += p_j.mass_ * (viscous_force_vec + pressure_force_vec);
        p_j.a_vec_ -= p_i.mass_ * (viscous_force_vec + pressure_force_vec);
        return;
    }
    return;
};

void internalForce(ParticleGroup& particle_group, const double h, const int dim) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        for (int j = 0; j < i+1; ++j) {
            internalForce(particle_group[i], particle_group[j], h, dim);
        }
    }
    return;
};

void wallForce(Particle& p_i, Particle& p_w, const double h, const int dim) {
    const Container::DoubleArray delta_x_iw_vec = p_i.x_vec_ - p_w.x_vec_;
    const double r_iw = Container::norm(delta_x_iw_vec);
    if (r_iw > radius_ratio * h) {
        return;
    }
    else {
        double coefficient = p_i.mass_ * p_i.p_;
        const Container::DoubleArray kernel_value_gradient = kernelValueGradient(delta_x_iw_vec, h, dim);
        if (p_i.p_ < 0) {
            return;
        }
        else if (p_i.rho_ > ConstVariable::Physics::Water::rho_0) {
            coefficient *= -2.0 / p_i.rho_ / p_i.rho_;
            p_i.a_vec_ += coefficient * kernel_value_gradient;
            return;
        }
        else {
            coefficient *= - 1/p_i.rho_/p_i.rho_ - 1/ConstVariable::Physics::Water::rho_02;
            p_i.a_vec_ += coefficient * kernel_value_gradient;
            return;
        }
    }
};

void wallForce(ParticleGroup& particle_group, ParticleGroup& wall_particle_group, const double h, const int dim) {
    if (particle_group.type_ == ParticleType::Wall) {
        return;
    }
    for (int i = 0; i < particle_group.n_particles_; ++i) {
        for (int j = 0; j < wall_particle_group.n_particles_; ++j) {
            wallForce(particle_group[i], wall_particle_group[j], h, dim);
        }
    }
    return;
};

};

#endif // INTERACTIVE_ACTION_CPP_