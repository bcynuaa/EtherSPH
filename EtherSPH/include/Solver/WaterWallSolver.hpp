/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-10 15:23:26
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: header file for WaterWallSolver
 */

#ifndef WATER_WALL_SOLVER_CPP_
#define WATER_WALL_SOLVER_CPP_

#include <iostream>
#include <string>
#include <fstream>
#include "ConstVariable/ConstVariable.hpp"
#include "Container/Container.hpp"
#include "KernelFunction/KernelFunction.hpp"
#include "Particle/Particle.hpp"

namespace Solver {

struct WaterWallSolver {
    const int dim_;
    const double h_;
    const double dt_;
    const Container::DoubleArray body_force_vec_;
    const int total_step_;
    const int output_step_;
    Particle::ParticleGroup& water_particle_group_;
    Particle::ParticleGroup& wall_particle_group_;
    std::string output_path_;
    std::string output_file_name_;

    WaterWallSolver(
        const int dim,
        const double h,
        const double dt,
        Container::DoubleArray& body_force_vec,
        const int total_step,
        const int output_step,
        Particle::ParticleGroup& water_particle_group,
        Particle::ParticleGroup& wall_particle_group,
        std::string output_path,
        std::string output_file_name
    ) : dim_(dim),
        h_(h),
        dt_(dt),
        body_force_vec_(body_force_vec),
        total_step_(total_step),
        output_step_(output_step),
        water_particle_group_(water_particle_group),
        wall_particle_group_(wall_particle_group), 
        output_path_(output_path),
        output_file_name_(output_file_name){
            return;
        };

    std::string outputFileNameAtStep(const int i_step) {
        return output_path_ + output_file_name_ + std::to_string(i_step) + ".vtp";
    };

    void writeCurrentStepToVtp(const int i_step) {
        const std::string file_name = outputFileNameAtStep(i_step);
        std::ofstream file(file_name);
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "<PolyData>\n";
        file << "<Piece NumberOfPoints=\"" << water_particle_group_.n_particles_ + wall_particle_group_.n_particles_ << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
        file << "<Points>\n";
        file << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < water_particle_group_.n_particles_; ++i) {
            file << water_particle_group_[i].x_vec_[0] << " " << water_particle_group_[i].x_vec_[1] << " " << 0. << "\n";
        }
        for (int i = 0; i < wall_particle_group_.n_particles_; ++i) {
            file << wall_particle_group_[i].x_vec_[0] << " " << wall_particle_group_[i].x_vec_[1] << " " << 0. << "\n";
        }
        file << "</DataArray>\n";
        file << "</Points>\n";
        file << "<PointData>\n";
        file << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < water_particle_group_.n_particles_; ++i) {
            file << water_particle_group_[i].v_vec_[0] << " " << water_particle_group_[i].v_vec_[1] << " " << 0. << "\n";
        }
        for (int i = 0; i < wall_particle_group_.n_particles_; ++i) {
            file << wall_particle_group_[i].v_vec_[0] << " " << wall_particle_group_[i].v_vec_[1] << " " << 0. << "\n";
        }
        file << "</DataArray>\n";
        file << "<DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (int i = 0; i < water_particle_group_.n_particles_; ++i) {
            file << water_particle_group_[i].p_ << "\n";
        }
        for (int i = 0; i < wall_particle_group_.n_particles_; ++i) {
            file << wall_particle_group_[i].p_ << "\n";
        }
        file << "</DataArray>\n";
        file << "<DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (int i = 0; i < water_particle_group_.n_particles_; ++i) {
            file << water_particle_group_[i].rho_ << "\n";
        }
        for (int i = 0; i < wall_particle_group_.n_particles_; ++i) {
            file << wall_particle_group_[i].rho_ << "\n";
        }
        file << "</DataArray>\n";
        file << "</PointData>\n";
        file << "</Piece>\n";
        file << "</PolyData>\n";
        file << "</VTKFile>\n";
        file.close();
    };

    void solve() {
        if (output_path_.back() != '/') {
            output_path_ += '/';
        }
        std::string command = "mkdir -p " + output_path_;
        system(command.c_str());
        // solve
        int count_write_step = 0;
        for (int i_step = 0; i_step < this->total_step_+1; ++i_step) {
            if (i_step % output_step_ == 0) {
                std::cout << "step: " << i_step << std::endl;
                writeCurrentStepToVtp(count_write_step);
                ++count_write_step;
            }
            Particle::balanceMass(this->water_particle_group_, this->h_, this->dim_);
            Particle::updateDensity(this->water_particle_group_, this->dt_);
            Particle::updatePressure(this->water_particle_group_, this->dt_);
            Particle::internalForce(this->water_particle_group_, this->h_, this->dim_);
            Particle::wallForce(this->water_particle_group_, this->wall_particle_group_, this->h_, this->dim_);
            Particle::updateVelocity(this->water_particle_group_, this->dt_, this->body_force_vec_);
            Particle::updatePosition(this->water_particle_group_, this->dt_);
        }
    };
};

};

#endif // WATER_WALL_SOLVER_CPP_