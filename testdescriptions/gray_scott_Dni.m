function [Dni] = gray_scott_Dni(t,dt,U,V)
  % Usage: [Dni] = gray_scott_Dni(t,dt,U,V)
  %
  % Gray Scott test problem.  The state vector is U = [u,v]^T.  The problem is:
  %   udot = Du \nabla^2 u + A*(1 - u) - u*v*v
  %   vdot = Dv \nabla^2 v + u*v*v - (A + B)*v
  % <=>
  %   U' = F(U)
  % on a periodic rectangular domain (x,y) \in \Omega.
  %
  % The values of the parameters A, B, Du and Dv, as well as the
  % space-time interval \Omega \times (t0,tf], and initial condition
  % U(t0,x,y) are all defined in the file gray_scott_setup.m
  %
  % This routine implements the quantity
  %   Dni = F(V) - Jac*V - Nn
  % where
  %   Jac  = F_U(U)
  %   Nn   = F(U) - F_U(U)*U
  %
  % Here,
  %   F(V) = [ Du \nabla^2 V_u + A - A*V_u - V_u*V_v*V_v     ]
  %          [ Dv \nabla^2 V_v + V_u*V_v*V_v - A*V_v - B*V_v ]
  %   F_U(U) = [ Du \nabla^2 - A - U_v*U_v,  -2*U_u*U_v                     ]
  %            [    U_v*U_v,                 Dv \nabla^2 - (A+B) + 2*U_u*U_v]
  %   F_U(U)V = [ Du \nabla^2 V_u - A*V_u - U_v*U_v*V_u -2*U_u*U_v*V_v]                     ]
%               [ U_v*U_v*V_u + Dv \nabla^2 V_v - (A+B)*V_v + 2*U_u*U_v*V_v]
  %   Nn(U) = [A + 2*U_u*U_v*U_v]
  %           [  -2*U_u*U_v*U_v ]
  % so combining these, we have
  %   Dni = [ -V_u*V_v*V_v + U_v*U_v*V_u + 2*U_u*U_v*V_v - 2*U_u*U_v*U_v ]
  %         [  V_u*V_v*V_v - U_v*U_v*V_u - 2*U_u*U_v*V_v + 2*U_u*U_v*U_v ]
  %
  % Daniel R. Reynolds
  % February 15, 2022
  %
  % adapted from code by John Loffeld

  % access problem parameters via global variables
  global nx
  global ny
  global idx

  % initialize Dni
  Dni = zeros(size(U));

  % compute Dni
  for j=1:ny
    for i=1:nx
      U_u = U(idx(1,i,j));
      U_v = U(idx(2,i,j));
      V_u = V(idx(1,i,j));
      V_v = V(idx(2,i,j));
      Dni(idx(1,i,j)) = -V_u*V_v*V_v + U_v*U_v*V_u + 2*U_u*U_v*V_v - 2*U_u*U_v*U_v;
      Dni(idx(2,i,j)) = V_u*V_v*V_v - U_v*U_v*V_u - 2*U_u*U_v*V_v + 2*U_u*U_v*U_v;
    end
  end

end
