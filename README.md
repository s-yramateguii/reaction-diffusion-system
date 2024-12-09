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
- $$D_1=0.1=D_2$$ are the diffusion coefficients for $$u$$ and $$v$$, respectively.
- $$\nabla^2 = \partial_x^2 + \partial_y^2$$

The equations describe the evolution of the concentration of two species (or fields), with the terms for diffusion, reaction, and nonlinearity. The simulation uses both spectral methods (periodic boundary conditions) for solving the system of equations in Fourier space and Chebyshev polynomials (no-flux boundary conditions) in real space.

## Initial Conditions

### Spatial Domain:
The simulation is carried out on a 2D grid with the spatial domain defined as:
- $$L_x = 20$$ (width of the domain)
- $$L_y = 20$$ (height of the domain)
- The grid is discretized into $$64 \times 64$$ points

### Field Initialization:
The initial fields $$u(x, y)$$ and $$v(x, y)$$ are set to represent a spiral pattern.

$$
u(x, y) = \tanh(r) \cdot \cos(m\theta - r)
$$

$$
v(x, y) = \tanh(r) \cdot \sin(m\theta - r)
$$

Where:
- $$r = \sqrt{x^2 + y^2}$$ is the radial distance.
- $$\theta = \text{angle}(x + iy)$$ is the angle in polar coordinates.
- $$m$$ is the winding number that determines the number of arms of the spiral.

## Project Components

### 1. Fourier Spectral Method (Part A)
In the first part, the system is solved using the Fourier spectral method, where the solution is transformed to Fourier space using the Fast Fourier Transform (FFT), allowing for efficient handling of linear diffusion and reaction terms. The nonlinear terms are handled directly in real space, while the linear diffusion terms are treated in Fourier space.

### 2. Chebyshev Spectral Method (Part B)
In the second part, a Chebyshev spectral method is used, involving Chebyshev polynomials for spatial discretization. The diffusion operators are approximated using Chebyshev differentiation matrices, and boundary conditions are applied using no-flux boundary conditions.

### 3. Animation
Both methods generate time-series solutions for $$U$$ and $$V$$, and these solutions are visualized in real-time using animations.
