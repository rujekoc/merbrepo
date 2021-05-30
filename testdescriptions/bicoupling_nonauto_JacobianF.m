function [J] = bicoupling_nonauto_JacobianF(t,U)
  % J = bicoupling_nonauto_JacobianF(t,U)
  % Jacobian for fast portion (fixed linearization)
  % of non-autonomous bidirectional coupling test problem
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

  n = length(U);
  J = zeros(n,n);
  x = U(1); y = U(2); z = U(3);
  J(1,2) = w;
  J(2,1) = -w;

end
