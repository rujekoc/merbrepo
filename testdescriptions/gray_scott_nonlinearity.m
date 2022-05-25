function [Udot] = gray_scott_nonlinearity(t,U)
  % Udot = gray_scott_nonlinearity(t,U)
  % Nonlinearity for Gray Scott test problem. The state vector is U = [u,v]^T.  The problem is:
  %   udot = Du \nabla^2 u + A*(1 - u) - u*v*v
  %   vdot = Dv \nabla^2 v + u*v*v - (A + B)*v
  % <=>
  %   udot = Du \nabla^2 u + A - Au - u*v*v
  %   vdot = Dv \nabla^2 v + u*v*v - Av - Bv
  % <=>
  %   U' = F(U)
  % on a periodic rectangular domain (x,y) \in \Omega.
  %
  % The values of the parameters A, B, Du and Dv, as well as the
  % space-time interval \Omega \times (t0,tf], and initial condition
  % U(t0,x,y) are all defined in the file gray_scott_setup.m
  %
  % This routine implements the quantity
  %   Nn = F(U) - F_U(U)*U
  % Here, since
  %   F_U(U) = [ Du \nabla^2 - A - v*v,  -2*u*v                     ]
  %            [    v*v,                 Dv \nabla^2 - (A+B) + 2*u*v]
  % then
  %   F_U(U)*U = [ Du \nabla^2 u - Au - 3*u*v*v    ]
  %              [ Dv \nabla^2 v - (A+B)v + 3*u*v*v]
  % and so
  %   Nn = [A + 2*u*v*v]
  %        [-2*u*v*v]
  % Daniel R. Reynolds
  % February 15, 2022
  %
  % adapted from code by John Loffeld

  % access problem parameters via global variables
  global nx
  global ny
  global A
  global idx

  % initialize Udot
  Udot = zeros(size(U));

  % compute residual nonlinearity
  for j=1:ny
    for i=1:nx
      u = U(idx(1,i,j));
      v = U(idx(2,i,j));
      Udot(idx(1,i,j)) = A + 2*u*v*v;
      Udot(idx(2,i,j)) = -2*u*v*v;
    end
  end

end
