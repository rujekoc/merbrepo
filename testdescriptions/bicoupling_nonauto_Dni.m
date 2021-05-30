function [Dni] = bicoupling_nonauto_Dni(t,dt,U,V)
  % Udot = bicoupling_nonauto_Dni(t,dt,U,V)
  % non-autonomous bidirectional coupling test problem
  %
  %   [x]'   [wy-z-pt]
  %   [y]  = [-wx]
  %   [z]    [-lz - lpt - p(x - (az)/(al+bw) - (apt)/(al + bw))^2 ...
  %                     - p(y - (bz)/(al+bw) - (bt)/(al + bpw))^2]
  % <=>
  %   U' = F(t,U)
  %
  % This routine implements the quantity
  %   Dni = F(t+dt,V) - Jac*V - Nn - dt*JacT
  % where
  %   Jac  = F_U(t,U)
  %   JacT = F_t(t,U)
  %   Nn   = F(t,U) - F_U(t,U)*U

  % Daniel R. Reynolds
  % June 2020

  global a
  global b
  global w
  global l
  global p

  n = length(U);
  Dni = zeros(n,1);
  Ux = U(1); Uy = U(2); Uz = U(3);
  Vx = V(1); Vy = V(2); Vz = V(3);
  tau = t+dt;
  xfac = (Ux - a*Uz/(a*l + b*w) - a*p*t/(a*l + b*w));
  yfac = (Uy - b*Uz/(a*l + b*w) - b*p*t/(a*l + b*w));


  Dni(3) = p*( - (Vx - a*Vz/(a*l + b*w) - a*p*tau/(a*l + b*w))^2 ...
               - (Vy - b*Vz/(a*l + b*w) - b*p*tau/(a*l + b*w))^2 ...
               + xfac*xfac + yfac*yfac ...
               + 2*xfac*(Vx-Ux) ...
               + 2*yfac*(Vy-Uy) ...
               + (2*a)/(a*l + b*w)*xfac*(Uz-Vz) ...
               + (2*b)/(a*l + b*w)*yfac*(Uz-Vz) ...
               - dt*(2*a*p)/(a*l + b*w)*xfac  ...
               - dt*(2*b*p)/(a*l + b*w)*yfac );

end
