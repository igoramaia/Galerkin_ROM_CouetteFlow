# Galerkin ROMs for Couette flow
Galerkin_ROM_Couette_flow is a set of codes to generate, run and explore Galerkin-based reduced order models (ROMs) for Couette flow. Such ROMs can be used for multiple purposes, such as estimating flow turbulence statistics [1], deriving generalised-quasilinear approximations [2] and performing control optimisation [3]. The models use controllability modes as basis functions, which are derived from the reponse of the linearised Navier-Stokes operator to white-noise forcing, as laid out in [1].

# The codes <h3>
ROM generation, time integration and post-processing are performed in the Matlab environment. The codes use a Fourier-Chebyshev discretisation of the domain, and a standard 4th/5th Runge-Kutta method for the time integration. They can be downloaded from the terminal by typing:

> `git clone https://github.com/igoramaia/Galerkin_ROM_CouetteFlow.git`

The files in the root are the main scripts to perform different tasks, which are described in the following. Auxiliary functions to perform mathematical operations are located in the **aux_functions** folder. The main scripts execute the following tasks:

* genROM.m: generates the modal basis and pre-computes the ROM operators, which are then stored in a .mat file. The main user parameters that can be altered are: the number of non-zero streamwise and spanwise wavenumbers included in the basis, the domain size and grid spacing. The latter should be carefully adapted in order to avoid alising issues.
* runROM.m: loads a given generated ROM and performs the time integration of the system of ODE's for a given Reynolds number and initial condition defined by the user. The time history of the coefficients and the modal basis are then used to reconstruct the full velocity field, $u(x,y,z,t)$, and to compute first- and second-order flow statistics.
* GQL_ROM.m: performs generalised-quasilinear (GQL) approximations in the ROM framework, based on a controllability criterion, following the mathematical formulation of [2]. The user has the choice of performing a standard GQL approximation, or a driven-GQL approximation. The linearisation procedure, and the corresponding nonlinear term in the ROM, are adapted according to this choice. Flow statistics and errors with respect to the original ROM are computed in the end.
* forc_optim_ROM.m: this codes considers a ROM modified with the addition of a body forcing term composed of Stokes modes. The modal coefficients of the forcing are optmised to reduce the total fluctuation energy of the flow, following a gradient-descent algorithm, as laid out in [3].
* postproc_optim_forc.m: loads the optimal forcing coefficients, reconstructs the optimal forcing and performs simulations of the forced flow. Different aspects of the forced flow are explored, such as its energy budget, linear stability and transient growth.

# References <h3>
The codes are based on the mathematical formulation developed in the following articles. Please cite them accordingly:

1. A. V. G. Cavalieri and P. A. S. Nogueira. "Reduced-order Galerkin models of plane Couette flow". Physical Reviw Fluids, 7, L102601, 2022. https://doi.org/10.1103/PhysRevFluids.7.L102601
2. I. A. Maia and A. V. G. Cavalieri. "Modal-based generalised quasilinear approximations for turbulent plane Couette flow". Theoretical and Computational Fluid Dynamics, 38, 313-330, 2024. https://doi.org/10.1007/s00162-024-00691-4
3. I. A. Maia and A. V. G. Cavalieri. "Turbulence suppression in plane Couette flow using reduced-order models". Journal of Fluid Mechanics, 1014, A18, 2025. https://doi.org/10.1017/jfm.2025.10258
