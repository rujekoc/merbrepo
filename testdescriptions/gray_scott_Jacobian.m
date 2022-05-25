function [J] = gray_scott_Jacobian(t,U)
  % Usage: [J] = gray_scott_Jacobian(t,U)
  %
  % Gray Scott test problem RHS Jacobian routine
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
  global xl
  global xr
  global yl
  global yr
  global dx
  global dy
  global A
  global B
  global Du
  global Dv
  global idx
  global neighbors

  % initialize J
  J = spalloc(2*nx*ny,2*nx*ny,12*nx*ny);

  % compute diffusion contributions
  uxcoef = Du/(dx*dx);
  uycoef = Du/(dy*dy);
  vxcoef = Dv/(dx*dx);
  vycoef = Dv/(dy*dy);
  for j=1:ny
    for i=1:nx
      nb = neighbors(i,j);
      il = nb(1); ir = nb(2); jl = nb(3); jr = nb(4);
      J(idx(1,i,j),idx(1,ir,j)) = J(idx(1,i,j),idx(1,ir,j)) + uxcoef;
      J(idx(1,i,j),idx(1,il,j)) = J(idx(1,i,j),idx(1,il,j)) + uxcoef;
      J(idx(1,i,j),idx(1,i,jr)) = J(idx(1,i,j),idx(1,i,jr)) + uycoef;
      J(idx(1,i,j),idx(1,i,jl)) = J(idx(1,i,j),idx(1,i,jl)) + uycoef;
      J(idx(1,i,j),idx(1,i,j))  = J(idx(1,i,j),idx(1,i,j))  - 2*(uxcoef+uycoef);

      J(idx(2,i,j),idx(2,ir,j)) = J(idx(2,i,j),idx(2,ir,j)) + vxcoef;
      J(idx(2,i,j),idx(2,il,j)) = J(idx(2,i,j),idx(2,il,j)) + vxcoef;
      J(idx(2,i,j),idx(2,i,jr)) = J(idx(2,i,j),idx(2,i,jr)) + vycoef;
      J(idx(2,i,j),idx(2,i,jl)) = J(idx(2,i,j),idx(2,i,jl)) + vycoef;
      J(idx(2,i,j),idx(2,i,j))  = J(idx(2,i,j),idx(2,i,j))  - 2*(vxcoef+vycoef);
    end
  end

  % compute reaction contributions
  for j=1:ny
    for i=1:nx
      u = U(idx(1,i,j));
      v = U(idx(2,i,j));
      J(idx(1,i,j),idx(1,i,j)) = J(idx(1,i,j),idx(1,i,j)) - A - v*v;
      J(idx(1,i,j),idx(2,i,j)) = J(idx(1,i,j),idx(2,i,j)) - 2*u*v;

      J(idx(2,i,j),idx(1,i,j)) = J(idx(2,i,j),idx(1,i,j)) + v*v;
      J(idx(2,i,j),idx(2,i,j)) = J(idx(2,i,j),idx(2,i,j)) + 2*u*v - (A + B);
    end
  end

end
