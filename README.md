# MERB: Multirate Exponential Rosenbrock methods #

#### Rujeko Chinomona<sup>1</sup>, Daniel R. Reynolds<sup>1</sup>, and Vu Thai Luan<sup>2</sup> ####
<sup>1</sup>Department of Mathematics, Southern Methodist University |
<sup>2</sup>Department of Mathematics, Mississippi State University

This repository contains MATLAB files for test problems highlighted in the paper, [Vu Thai Luan, Rujeko Chinomona & Daniel R. Reynolds, "Multirate exponential Rosenbrock methods," arXiv:2106.05385, 2021](http://arxiv.org/abs/2106.05385).

Multirate Exponential Rosenbrock (MERB) methods of orders three (`MERB3`), four (`MERB4`), five (`MERB5`),and six (MERB6) are implemented on the additively split multirate problem:

  u' = F(t,u) = J<sub>n</sub>u + V<sub>n</sub>t + N<sub>n</sub>(t,u)

where J<sub>n</sub> is the Jacobian of F at (t<sub>n</sub>,u<sub>n</sub>), V<sub>n</sub> is the partial derivative of F with respect to time at (t<sub>n</sub>,u<sub>n</sub>), and N<sub>n</sub>(t,u) = F(t,u) - J<sub>n</sub>u - V<sub>n</sub>t.

MERB methods by design require an update of the Jacobian at each time step, an approach we refer to as dynamic linearization. The fast linear part of the algorithm is F<sub>f</sub>(t,u) = J<sub>n</sub>u and the slow
part is F<sub>s</sub>(t,u) = V<sub>n</sub>t + N<sub>n</sub>(t,u).

We run comparison tests with Multirate Exponential Runge-Kutta (MERK) methods ( third-order ```MERK3```, fourth-order ```MERK4```, and fifth-order ```MERK5```) from [V.T. Luan, R. Chinomona and D.R. Reynolds, SIAM Journal on Scientific Computing, 2020](https://doi.org/10.1137/19M125621X) and Multirate Infinitesimal Generalized-structure Additive Runge-Kutta (MRIGARK) methods (third-order ```MRI-GARK-ERK33a``` and fourth-order ```MRI-GARK-ERK45a```) from [A. Sandu, SIAM Journal on Numerical Analysis, 2019](https://doi.org/10.1137/18M1205492).

MERK and MRIGARK methods use both dynamic linearization and fixed linearization fast-slow splittings. For each test problem optimal time scale separation factors m, which are integer ratios between the slow and fast time step sizes are determined and used to compare methods at their peak efficiency.

The two test problems included in this repository are a 1D reaction-diffusion problem and a non-autonomous bidirectional coupling problem with descriptions in `testdescriptions`. An analytical solution for the bidirectional coupling problem is provided, while a reference solution for the reaction-diffusion problem is computed by ```ode15s```. The drivers for the considered methods are included in the folders `merb`, `merk`, `mrigark`. For each method and each test problem the drivers output maximum absolute errors, root-mean-square errors, rate of convergence results, number of slow and fast function evaluations, and runtimes.
