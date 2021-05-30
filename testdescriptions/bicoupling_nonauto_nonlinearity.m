function [Udot] = bicoupling_nonauto_nonlinearity(t,U)
  % Udot = bicoupling_nonauto_nonlinearity(t,U)
  % Nonlinearity for non-autonomous bidirectional coupling test problem
  %
  %   [x]'   [wy-z-pt]
  %   [y]  = [-wx]
  %   [z]    [-lz - lpt - p(x - (az)/(al+bw) - (apt)/(al + bw))^2 ...
  %                     - p(y - (bz)/(al+bw) - (bt)/(al + bpw))^2]
  % <=>
  %   [x]'    [x]
  %   [y]  = J[y] + N(U,t)
  %   [z]     [z]
  % Where J is the Jacobian of the right-hand side with respect to U.
  % This routine implements N(U,t).
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

  n = length(U);
  Udot = zeros(n,1);
  x = U(1); y = U(2); z = U(3);

  xfac = (x - a*z/(a*l + b*w) - a*p*t/(a*l + b*w));
  yfac = (y - b*z/(a*l + b*w) - b*p*t/(a*l + b*w));
  Udot(1) = -p*t;
  Udot(2) = 0;
  Udot(3) = p*(-l*t + xfac*(2*x - xfac - (2*a)/(a*l + b*w)*z) ...
                    + yfac*(2*y - yfac - (2*b)/(a*l + b*w)*z));
end
