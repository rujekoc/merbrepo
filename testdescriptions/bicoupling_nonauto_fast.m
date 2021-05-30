function [Udot] = bicoupling_nonauto_fast(t,U)
  % Udot = bicoupling_nonauto_fast(t,U)
  % (fixed linearization) Computes fast rhs for
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
  Udot(1) = w*y;
  Udot(2) = -w*x;

end
