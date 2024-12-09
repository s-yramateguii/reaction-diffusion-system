# Reaction-Diffusion System

This repository simulates a two-component reaction-diffusion system, commonly used to model pattern formation in biological and chemical systems. The simulation employs both Fourier and Chebyshev polynomials to solve the system of partial differential equations.

## Overview

The system consists of two fields, $$U$$ and $$V$$, which interact according to the following nonlinear partial differential equations:

$$
\frac{\partial U}{\partial t} = D_1 \nabla^2 U + \lambda(U, V) U - \omega(U, V) V
$$

$$
\frac{\partial V}{\partial t} = D_2 \nabla^2 V + \lambda(U, V) V + \omega(U, V) U
$$

Where:
- $$\lambda(U, V) = 1 - (U^2 + V^2)$$
- $$\omega(U, V) = -\beta (U^2 + V^2)$$

The equations describe the evolution of the concentration of two species (or fields), with the terms for diffusion, reaction, and nonlinearity. The simulation uses both spectral methods (periodic boundary conditions) for solving the system of equations in Fourier space and Chebyshev polynomials (no-flux boundary conditions) in real space.

## Project Components

### 1. Fourier Spectral Method (Part A)
In the first part, the system is solved using the Fourier spectral method, where the solution is transformed to Fourier space using the Fast Fourier Transform (FFT), allowing for efficient handling of linear diffusion and reaction terms. The nonlinear terms are handled directly in real space, while the linear diffusion terms are treated in Fourier space.

### 2. Chebyshev Spectral Method (Part B)
In the second part, a Chebyshev spectral method is used, involving Chebyshev polynomials for spatial discretization. The diffusion operators are approximated using Chebyshev differentiation matrices, and boundary conditions are applied using no-flux boundary conditions.

### 3. Animation
Both methods generate time-series solutions for $$U$$ and $$V$$, and these solutions are visualized in real-time using animations.
