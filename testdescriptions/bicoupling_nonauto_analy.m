function U = bicoupling_nonauto_analy(t)
  % U = bicoupling_nonauto_analy(t)
  % Analytical solution of non-autonomous bidirectional coupling test problem
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

  wt = vpa(w*t,20);
  elt = exp(vpa(-l*t,20));
  U = double([cos(wt) + a*elt; -sin(wt) + b*elt; (a*l + b*w)*elt - p*t]);

end
