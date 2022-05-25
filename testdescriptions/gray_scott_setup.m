function [U0,t0,tf] = gray_scott_setup()
  % Usage: [U0,t0,tf] = gray_scott_setup()
  %
  % Gray-Scott test problem setup routine.  This routine
  % generates the initial condition U(t0) = U0 for the Gray-Scott
  % test problem, as well as the time interval over which it should
  % be solved, [t0,tf].  Additionally, this routine defines some
  % test problem parameters as global variables; these are used by
  % the RHS routine(s) for this problem.  The problem is formulated
  % as follows.
  %
  % The state vector is U = [u,v]^T.  The problem is:
  %   udot = Du \nabla^2 u - u*v*v + A*(1 - u)
  %   vdot = Dv \nabla^2 v + u*v*v - (A + B)*v
  % on the domain (x,y) \in [0,L]^2, t \in (0,0.2], with parameters
  %   A = 0.625
  %   B = 0.25
  %   Du = 0.0312
  %   Dv = 0.0156
  % The initial condition U(t0,x,y) is defined as:
  %   u(t0,x,y) = 1 - exp(-150*((x-L/2)^2 + (y-L/2)^2))
  %   v(t0,x,y) = exp(-150*((x-L/2)^2 + 2*(y-L/2)^2))
  % and the boundary conditions are periodic in both directions.


  % define reusable global variables
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

  % set problem parameters
  nx = 50;
  ny = 50;
  xl = 0;
  xr = 1;
  yl = 0;
  yr = 1;
  dx = (xr-xl)/nx;
  dy = (yr-yl)/ny;
  A  = 0.625;
  B  = 0.25;
  Du = 0.312;
  Dv = 0.156;

  % set time interval
  t0 = 0;
  tf = 0.2;

  % 2D -> 1D mapping utility function for physical -> logical indexing
  idx = @(l,i,j) (((j-1)*ny + (i-1))*2 + l);

  % neighbor index mapping function to handle periodicity
  neighbors = @(i,j) [(i>1)*(i-1)+(i==1)*nx,(i<nx)*(i+1)+(i==nx)*1,(j>1)*(j-1)+(j==1)*ny,(j<ny)*(j+1)+(j==ny)*1];

  % set initial condition
  U0 = zeros(2*nx*ny,1);
  for j=1:ny
    y = yl + (j-0.5)*dy;
    for i=1:nx
      x = xl + (i-0.5)*dx;
      U0(idx(1,i,j)) = 1 - exp(-150*((x-xr/2)^2 + (y-yr/2)^2));
      U0(idx(2,i,j)) = exp(-150*((x-xr/2)^2 + 2*(y-yr/2)^2));
    end
  end

end
