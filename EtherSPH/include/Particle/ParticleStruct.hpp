/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-07 21:13:08
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: header file for ParticleStruct
 */

#ifndef PARTICLE_STRUCT_HPP_
#define PARTICLE_STRUCT_HPP_

#include <iostream>
#include "Container/Container.hpp"
#include "KernelFunction/KernelFunction.hpp"

namespace Particle {

enum ParticleType {
    Water= 0,
    Wall = 1
};

struct Particle {
    Container::DoubleArray x_vec_;
    Container::DoubleArray v_vec_;
    Container::DoubleArray a_vec_;
    double mass_;
    double rho_;
    double p_;
    double drho_;

    Particle(const int dim) {
        this->x_vec_ = Container::DoubleArray(dim);
        this->v_vec_ = Container::DoubleArray(dim);
        this->a_vec_ = Container::DoubleArray(dim);
        this->mass_ = 0.0;
        this->rho_ = ConstVariable::Physics::Water::rho_0;
        this->p_ = 0.0;
        this->drho_ = 0.0;
    };

    Particle() : Particle(2) {};

    Particle(const int dim, const ParticleType type) {
        this->x_vec_ = Container::DoubleArray(dim);
        this->v_vec_ = Container::DoubleArray(dim);
        this->a_vec_ = Container::DoubleArray(dim);
        this->mass_ = 0.0;
        this->rho_ = ConstVariable::Physics::Water::rho_0;
        this->p_ = 0.0;
        this->drho_ = 0.0;
    };

    Particle(const Particle& particle) {
        Container::allocate(this->x_vec_, particle.x_vec_);
        Container::allocate(this->v_vec_, particle.v_vec_);
        Container::allocate(this->a_vec_, particle.a_vec_);
        this->mass_ = particle.mass_;
        this->rho_ = particle.rho_;
        this->p_ = particle.p_;
        this->drho_ = particle.drho_;
    };

    ~Particle() {};

    Particle operator=(const Particle& particle) {
        Container::allocate(this->x_vec_, particle.x_vec_);
        Container::allocate(this->v_vec_, particle.v_vec_);
        Container::allocate(this->a_vec_, particle.a_vec_);
        this->mass_ = particle.mass_;
        this->rho_ = particle.rho_;
        this->p_ = particle.p_;
        this->drho_ = particle.drho_;
        return *this;
    };

    // overload operator <<
    friend std::ostream& operator<<(std::ostream& os, const Particle& particle) {
        os << "x_vec_: " << particle.x_vec_ << std::endl;
        os << "v_vec_: " << particle.v_vec_ << std::endl;
        os << "a_vec_: " << particle.a_vec_ << std::endl;
        os << "mass_: " << particle.mass_ << std::endl;
        os << "rho_: " << particle.rho_ << std::endl;
        os << "p_: " << particle.p_ << std::endl;
        os << "drho_: " << particle.drho_;
        return os;
    };
};

struct ParticleGroup {
    int dim_;
    int n_particles_;
    ParticleType type_;
    Container::Array<Particle> particles_;

    ParticleGroup(const int dim, const int n_particles, const ParticleType type) {
        this->dim_ = dim;
        this->n_particles_ = n_particles;
        this->type_ = type;
        this->particles_ = Container::Array<Particle>(n_particles);
        for (int i = 0; i < n_particles; ++i) {
            this->particles_[i] = Particle(dim);
        }
    };

    ParticleGroup(const ParticleGroup& particle_group) {
        this->dim_ = particle_group.dim_;
        this->n_particles_ = particle_group.n_particles_;
        this->type_ = particle_group.type_;
        this->particles_ = Container::Array<Particle>(particle_group.n_particles_);
        for (int i = 0; i < particle_group.n_particles_; ++i) {
            this->particles_[i] = particle_group.particles_[i];
        }
    };

    ~ParticleGroup() {};

    ParticleGroup operator=(const ParticleGroup& particle_group) {
        this->dim_ = particle_group.dim_;
        this->n_particles_ = particle_group.n_particles_;
        this->type_ = particle_group.type_;
        this->particles_ = Container::Array<Particle>(particle_group.n_particles_);
        for (int i = 0; i < particle_group.n_particles_; ++i) {
            this->particles_[i] = particle_group.particles_[i];
        }
        return *this;
    };

    // overload operator []
    Particle& operator[](const int index) {
        return this->particles_[index];
    };

    // overload operator <<
    friend std::ostream& operator<<(std::ostream& os, const ParticleGroup& particle_group) {
        os << "dim_: " << particle_group.dim_ << std::endl;
        os << "n_particles_: " << particle_group.n_particles_ << std::endl;
        os << "type_: " << particle_group.type_ << std::endl;
        os << "particles_: " << std::endl;
        for (int i = 0; i < particle_group.n_particles_; ++i) {
            os << std::endl << "particle " << i << ":" << std::endl;
            os << particle_group.particles_[i] << std::endl;
            os << "-------------------------------------------";
        }
        return os;
    };
};

};

#endif // PARTICLE_STRUCT_HPP_