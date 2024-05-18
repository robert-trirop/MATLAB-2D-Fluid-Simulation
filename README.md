# MATLAB 2D Fluid Simulation

This repository contains the MATLAB code for a fluid simulation, originally coded in 2017 and featured in [this YouTube video](https://www.youtube.com/watch?v=cM47L5RddsM).

<img src="images/simulation.png" alt="Fluid Simulation" width="400"/>

## Overview

The code simulates fluid flow using a combination of the finite difference method and the Runge-Kutta method for particle advection. The main functionalities include:
- Initialization of velocity, pressure fields, and particle positions.
- Application of a divergence-free condition to the velocity field using the Jacobi method.
- Self-advection of the velocity fields using the Lagrangian-Eulerian approach.
- Particle advection using the Runge-Kutta 4th order method.
- Visualization of particle positions over time to show the fluid flow.

## License

This project is licensed under the MIT License - see the `fluid.m` file for details.

## How to Use

1. Clone the repository.
2. Open `fluid.m` in MATLAB.
3. Run the script to see the fluid simulation.
4. Close the figure window to stop the simulation.

## Acknowledgments

This code was put together using knowledge I accumulated from the internet and my studies. If you use this code for your own animations, please give credit to the original author.
