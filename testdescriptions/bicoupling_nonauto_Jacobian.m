function [J] = bicoupling_nonauto_Jacobian(t,U)
  % J = bicoupling_nonauto_Jacobian(t,U)
  % Jacobian for non-autonomous bidirectional coupling test problem
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

  % Rujeko Chinomona
  % May 2020

  global a
  global b
  global w
  global l
  global p

  n = length(U);
  J = zeros(n,n);
  x = U(1); y = U(2); z = U(3);

  J(1,2) = w;
  J(1,3) = -1;
  J(2,1) = -w;
  J(3,1) = -2*p*(x - a*z/(a*l + b*w) - a*p*t/(a*l + b*w));
  J(3,2) = -2*p*(y - b*z/(a*l + b*w) - b*p*t/(a*l + b*w));
  J(3,3) = -l + (2*a*p)/(a*l + b*w)*(x - a*z/(a*l + b*w) - a*p*t/(a*l + b*w)) + ...
                (2*b*p)/(a*l + b*w)*(y - b*z/(a*l + b*w) - b*p*t/(a*l + b*w));

end
