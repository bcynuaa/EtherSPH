/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-10 17:19:29
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: test for WaterWallSolver
 */

#ifndef WATER_WALL_SOLVER_TEST_CPP_
#define WATER_WALL_SOLVER_TEST_CPP_

#include <iostream>
#include "Container/Container.hpp"
#include "Particle/Particle.hpp"
#include "Solver/Solver.hpp"

int main(int argc, char const *argv[]) {
    Container::DoubleArray gravity_vec = {0., -9.8};
    const double dt = 1e-4;
    const double h = 0.3;
    const double dy = 0.05;
    const int dim = 2;
    const int n_particles = 1200;
    Particle::ParticleGroup pg(dim, n_particles, Particle::ParticleType::Water);
    for (int i = 0; i < 30; ++i) {
        for (int j = 0; j < 40; ++j) {
            pg[i * 40 + j].x_vec_ = {2 * i * dy, -2 * j * dy};
            pg[i * 40 + j].v_vec_ = {0., 0.};
            pg[i * 40 + j].a_vec_ = {0., 0.};
            pg[i * 40 + j].mass_ = 4 * dy * dy * ConstVariable::Physics::Water::rho_0;
            pg[i * 40 + j].p_ = 2*ConstVariable::Physics::Water::rho_0 * 
                                ConstVariable::Physics::gravity * j * dy + 
                                ConstVariable::Physics::Water::p_0;
            pg[i * 40 + j].rho_ = ConstVariable::Physics::Water::rho_0 +
                                    (pg[i * 40 + j].p_ - ConstVariable::Physics::Water::p_0)
                                        / ConstVariable::Physics::Water::c_02;
        }
    }
    Particle::ParticleGroup wall_pg(dim, 500, Particle::ParticleType::Wall);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 100; ++j) {
            wall_pg[i * 100 + j].x_vec_ = {(2 * j - 65) * dy, -2 * (i + 40) * dy};
            wall_pg[i * 100 + j].v_vec_ = {0., 0.};
            wall_pg[i * 100 + j].a_vec_ = {0., 0.};
            wall_pg[i * 100 + j].mass_ = 4 * dy * dy;
            wall_pg[i * 100 + j].p_ = ConstVariable::Physics::Water::rho_0 * 
                                ConstVariable::Physics::gravity * j * dy + 
                                ConstVariable::Physics::Water::p_0;
            wall_pg[i * 100 + j].rho_ = ConstVariable::Physics::Water::rho_0 +
                                    (wall_pg[i * 100 + j].p_ - ConstVariable::Physics::Water::p_0)
                                        / ConstVariable::Physics::Water::c_02;
            wall_pg[i * 100 + j].mass_ *= wall_pg[i * 100 + j].rho_;
        }
    }
    const int total_step = 10000;
    const int output_step = 100;
    Solver::WaterWallSolver solver(dim, h, dt, gravity_vec, total_step, output_step, pg, wall_pg, "./water_wall_test/", "test");
    solver.solve();
    return 0;
};

#endif // WATER_WALL_SOLVER_TEST_CPP_