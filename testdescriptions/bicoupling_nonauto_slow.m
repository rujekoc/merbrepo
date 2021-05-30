function [Udot] = bicoupling_nonauto_slow(t,U)
  % Udot = bicoupling_nonauto_slow(t,U)
  % (fixed linearization) Computes slow rhs for
  % Non-autonomous bidirectional coupling test problem
  %
  % x' = wy-z-pt;
  % y' = -wx;
  % z' = -lz - lpt - p[x - (az)/(al + bw) - (apt)/(al + bw)]^2 - ...
  %        p[y - (bz)/(al+bw) - (bt)/(al + bpw)]^2;
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

  x = U(1); y = U(2); z = U(3);
  Udot = zeros(size(U));
  Udot(1) = -z-p*t;
  Udot(3) = -l*z - l*p*t - p*(x - a*z/(a*l + b*w) - a*p*t/(a*l + b*w))^2 - ...
                           p*(y - b*z/(a*l + b*w) - b*p*t/(a*l + b*w))^2;

end
