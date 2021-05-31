# MERB: Multirate Exponential Rosenbrock methods #

#### Rujeko Chinomona<sup>1</sup>, Daniel R. Reynolds<sup>1</sup>, and Vu Thai Luan<sup>2</sup> ####
<sup>1</sup>Department of Mathematics, Southern Methodist University |
<sup>2</sup>Department of Mathematics, Mississippi State University

This repository contains MATLAB files for test problems highlighted in the paper, [Vu Thai Luan, Rujeko Chinomona & Daniel R. Reynolds, "Multirate exponential Rosenbrock methods," arXiv:1904.06474, 2019](https://arxiv.org/abs/1904.06474).

*to insert arxiv link 

Multirate Exponential Rosenbrock (MERB) methods of orders three (`MERB3`), four (`MERB4`), five (`MERB5`),and six (MERB6) are implemented on the additively split multirate problem:

  u' = F(t,u) = J<sub>n</sub>u + V<sub>n</sub>t + N<sub>n</sub>(t,u)

where J<sub>n</sub> is the Jacobian of F at (t<sub>n</sub>,u<sub>n</sub>), V<sub>n</sub> is the partial derivative of F with respect to time at (t<sub>n</sub>,u<sub>n</sub>), and N<sub>n</sub>(t,u) = F(t,u) - J<sub>n</sub>u - V<sub>n</sub>t. MERB methods by design require an update of the Jacobian at each time step, an approach we refer to as dynamic linearization. The fast linear part of the algorithm is F<sub>f</sub>(t,u) = J<sub>n</sub>u and the slow
part is F<sub>s</sub>(t,u) = V<sub>n</sub>t + N<sub>n</sub>(t,u).

 *Comments on provision of Dni and Nn functions. 

We run comparison tests with Multirate Exponential Runge-Kutta (MERK) methods ( third-order ```MERK3```, fourth-order ```MERK4```, fifth-order ```MERK5```) from [V.T. Luan, R. Chinomona and D.R. Reynolds, SIAM Journal on Scientific Computing, 2020](https://doi.org/10.1137/19M125621X) and Multirate Infinitesimal Generalized-structure Additive Runge-Kutta (MRIGARK) methods (third-order ```MRI-GARK-ERK33a``` and fourth-order ```MRI-GARK-ERK45a```) from [A. Sandu, SIAM Journal on Numerical Analysis, 2019](https://doi.org/10.1137/18M1205492).

MERK and MRIGARK use both dynamic linearization and fixed linearization fast-slow splittings. For each test problem optimal time scale separation factors m are determined and use to compare methods at their peak efficiency.

There is a reaction-diffusion test problem and a bidirectional coupling problem.

For all test problems convergence stats, runtime stats, and convergence stats are printed.
