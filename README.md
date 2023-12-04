# EtherSPH

Let's break down to particles.

# About `SPH`

SPH method is mesh-free Lagrangian method to solve fluid dynamics. It is widely used in astrophysics, cosmology, and fluid dynamics. SPH is a method to solve the fluid dynamics by using particles. With fluid mass represented by particles, the PDE equation is solved by application of a 'kernel-interpolation' method.

# Why `EtherSPH`

Ether was a concept of a medium that permeated the universe and allowed light to propagate. However, it was disproved by the Michelson-Morley experiment.

For SPH method describes the fluid with 'particle' concept, and the way it approximates the derivatives of physical quantities is quite genuis but lacks of sufficient theoretical basis. Therefore, I named this project as `EtherSPH` to express the fact that SPH method is a kind of 'ether' in the field of fluid dynamics.

# Folder Structure

1. `EtherSPH`: The first `cpp` version of `EtherSPH`. It is a single-thread version without neighbor-searching acceleration. Although this version is not efficient and has bugs, it gives me a glance at this kind of method.
2. `EtherSPH-JL`: `Julia` is a high-level language with high efficiency. With better package management system and debugging tools, `Julia` is a good choice for scientific computing. This folder contains my trial on `Julia` version of `EtherSPH`. It contains the neighbour-searching acceleration and multi-dispatching for different kinds of particles, equation model and boundary conditions. I will keep on my work in this folder.
3. `EtherSPH-Notes`: This folder contains my notes on SPH method including my report slides on office meeting.