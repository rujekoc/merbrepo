function [J] = bicoupling_nonauto_JacobianT(t,U)
  % J = bicoupling_nonauto_JacobianT(t,U)
  % Derivative w.r.t time for non-autonomous bidirectional coupling test problem
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
  J = zeros(size(U));

  J(1) = -p;
  J(3) = -l*p + (2*a*p*p)/(a*l + b*w)*(x - a*z/(a*l + b*w) - a*p*t/(a*l + b*w)) + ...
                (2*b*p*p)/(a*l + b*w)*(y - b*z/(a*l + b*w) - b*p*t/(a*l + b*w));

end
