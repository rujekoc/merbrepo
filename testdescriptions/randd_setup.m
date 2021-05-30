function [U0,t0,tf] = randd_setup()
  % Reaction diffusion test problem from Savcenco et al.
  % "A multirate time stepping strategy for stiff ordinary differential equations", BIT 47 (2007)
  %      u' = ep*u_xx + u^2(1-u)
  % where u_x(0,t) = u_x(5,t) and u(x,0) = (1 + exp(lambda*(x-1))^(-1),
  % with parameters lambda = 5*sqrt(2).
  % We evaluate over the time interval [0,3].
  %
  % Rujeko Chinomona
  % Department of Mathematics
  % Southern Methodist University
  % May 2021

  % Set problem parameters
  global a
  global b
  global nx
  global ep
  global gamma
  global A

  a      = 0;
  b      = 5;
  nx     = 100;                              % number of spatial grid points
  x      = linspace(a,b,nx+1);
  ep     = 1e-2;
  gamma  = 1e-1;
  U0     = Y0_fun(x,gamma,ep);                            % initial condition

  % Set time parameters
  t0      = 0;
  tf      = 5;

  A = Afn(a,b,nx,ep);

end
function [A] = Afn(a,b,n,ep)
  % Usage: [A] = Afn(a,b,n,ep)
  % INPUTS : n - number of points in discretization
  %          a - left end-point
  %          b - right end-point
  %         ep - epsilon value
  %
  % Rujeko Chinomona
  % Department of Mathematics
  % Southern Methodist University

  % set spatial mesh size
  h = (b-a)/n;

  % create matrix storage
  A = spalloc(n+1,n+1,3*n+3);

  % interior rows
  for i=2:n
    A(i,i-1) = 1/h/h;
    A(i,i) = -2/h/h;
    A(i,i+1) = 1/h/h;
  end

  % boundary rows
  A(1,1:2) = [-2 2]/h^2;
  A(n+1,n:n+1) = [2 -2]/h^2;

  A = ep*A;
end
function Y = Y0_fun(x,gamma,ep)
  % usage: Y = Y0_fun(x,ep)
  % gives initial value for reaction diffusion.
  %
  % Rujeko Chinomona
  % Department of Mathematics
  % Southern Methodist University

  lambda = .5*sqrt(2*gamma/ep);
  Y = (1+ exp(lambda*(x-1))).^(-1);
  Y = Y';

end
