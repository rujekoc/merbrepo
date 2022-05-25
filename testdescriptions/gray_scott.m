function [Udot] = gray_scott(t,U)
  % Usage: [Udot] = gray_scott(t,U)
  %
  % Gray Scott test problem RHS routine
  %
  % The state vector is U = [u,v]^T.  The problem is:
  %   udot = Du \nabla^2 u + A*(1 - u) - u*v*v
  %   vdot = Dv \nabla^2 v + u*v*v - (A + B)*v
  % on a periodic rectangular domain (x,y) \in \Omega.
  %
  % The values of the parameters A, B, Du and Dv, as well as the
  % space-time interval \Omega \times (t0,tf], and initial condition
  % U(t0,x,y) are all defined in the file gray_scott_setup.m
  %
  % Daniel R. Reynolds
  % October 8, 2019
  %
  % adapted from code by John Loffeld

  % access problem parameters via global variables
  global nx
  global ny
  global dx
  global dy
  global A
  global B
  global Du
  global Dv
  global idx
  global neighbors

  % initialize Udot
  Udot = zeros(size(U));

  % compute diffusion contributions
  uxcoef = Du/(dx*dx);
  uycoef = Du/(dy*dy);
  vxcoef = Dv/(dx*dx);
  vycoef = Dv/(dy*dy);
  for j=1:ny
    for i=1:nx
      nb = neighbors(i,j);
      il = nb(1); ir = nb(2); jl = nb(3); jr = nb(4);
      Udot(idx(1,i,j)) = Udot(idx(1,i,j)) + uxcoef*(U(idx(1,ir,j))-2*U(idx(1,i,j))+U(idx(1,il,j))) ...
      + uycoef*(U(idx(1,i,jr))-2*U(idx(1,i,j))+U(idx(1,i,jl)));
      Udot(idx(2,i,j)) = Udot(idx(2,i,j)) + vxcoef*(U(idx(2,ir,j))-2*U(idx(2,i,j))+U(idx(2,il,j))) ...
      + vycoef*(U(idx(2,i,jr))-2*U(idx(2,i,j))+U(idx(2,i,jl)));
    end
  end

  % compute reaction contributions
  for j=1:ny
    for i=1:nx
      u = U(idx(1,i,j));
      v = U(idx(2,i,j));
      Udot(idx(1,i,j)) = Udot(idx(1,i,j)) + A*(1 - u) - u*v*v;
      Udot(idx(2,i,j)) = Udot(idx(2,i,j)) + u*v*v - (A + B)*v;
    end
  end

end
