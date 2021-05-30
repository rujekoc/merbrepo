function [Nn] = bicoupling_nonauto_nonlinearity_gen(tn,Un)
  % Nn = bicoupling_nonauto_nonlinearity_gen(tn,Un)
  % Generates nonlinearity function for non-autonomous bidirectional coupling test problem
  % depending on Jacobian evaluation
  %
  %   [x]'   [wy-z-pt]
  %   [y]  = [-wx]
  %   [z]    [-lz - lpt - p(x - (az)/(al+bw) - (apt)/(al + bw))^2 ...
  %                     - p(y - (bz)/(al+bw) - (bt)/(al + bpw))^2]
  % <=>
  %   [x]'    [x]
  %   [y]  = J[y] + N(t,U)
  %   [z]     [z]
  % Where J is the Jacobian of the right-hand side with respect to U.
  % This routine generates an anonymous function N(t,U) = F(t,U) - Jn*U.
  %
  % Tunable parameters: a,b,l,w,p
  % with aw = bl
  %
  % Analytical solution
  % x(t) = cos(wt) + ae^(-lt);
  % y(t) = -sin(wt) + be^(-lt);
  % z(t) = (al + bw)e^(-lt) - pt;

  global a
  global b
  global w
  global l
  global p

  xn = Un(1); yn = Un(2); zn = Un(3);
  afac = a/(a*l + b*w);
  bfac = b/(a*l + b*w);

  xfac_n = (xn - a*zn/(a*l + b*w) - a*p*tn/(a*l + b*w));
  yfac_n = (yn - b*zn/(a*l + b*w) - b*p*tn/(a*l + b*w));

  Nn = @(t,U) [-p*t; 0; -p*l*t - p*(U(1) - afac*(U(3) + p*t))^2 - ...
              p*(U(2) - bfac*(U(3) + p*t))^2 + ...
              2*p*xfac_n*(U(1) - afac*U(3)) + 2*p*yfac_n*(U(2) - bfac*U(3))];
end
