function [U0,t0,tf] = bicoupling_nonauto_setup(vars)
  % [U0,t0,tf] = bicoupling_nonauto_setup(vars)
  % Input: nothing or vector of global variables [a,b,w,l,p]
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
  %
  % Rujeko Chinomona
  % Department of Mathematics
  % Southern Methodist University
  % May 2021

  % define reusable global variables
  global a
  global b
  global w
  global l
  global p

  if (nargin == 1)
    % set problem parameters
    a = vars(1);
    b = vars(2);
    w = vars(3);
    l = vars(4);
    p = vars(5);

    if abs(a*w-b*l) > eps
      error('vars(1)*vars(3) should be equal to vars(2)*vars(4)');
      return;
    end

  else
    % set default problem parameters
    a = 1;
    b = 20;
    w = 100;
    l = 5;
    p = 1e-2;
  end

  % set time interval
  t0 = 0;
  tf = 1;

  U0 = bicoupling_nonauto_analy(0);

end
