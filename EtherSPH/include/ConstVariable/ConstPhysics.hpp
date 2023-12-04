/**
 * @ author: Chenyu Bao | bcynuaa@163.com
 * @ date: 2023-11-08 15:30:11
 * @ license: MIT
 * @ cxx standard: 11
 * @ description: header file for ConstPhysics
 */

#ifndef CONST_PHYSICS_HPP_
#define CONST_PHYSICS_HPP_

namespace ConstVariable {

namespace Physics {

const double gravity = 9.8;
const double avoid_singularity = 1e-2;

namespace Water {

const double rho_0 = 1000.;
const double rho_02 = rho_0 * rho_0;
const double c_0 = 50.;
const double c_02 = c_0 * c_0;
// const double p_0 = 1.013e5;
const double p_0 = 0.; // for convenience
const double mu_0 = 1e-3;
const double nu_0 = 1e-6;
// artificial compressibility to avoid accetional pressure oscillation
// Fluid mechanics and the SPH method: theory and applications, Violeau, Damien, 2012, p. 32
// Page 488
// however, this method is not used in this project
// i reference `SmoothedPartices.jl` to set the artificial compressibility
const double epsilon_rho = 1e-6;

};

};

};

#endif // CONST_PHYSICS_HPP_