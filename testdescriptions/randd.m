function [Udot] = randd(t,U)
  % Computes full right hand side of reaction diffusion problem

  global gamma
  global A
  Udot = A*U + gamma*(U.^2 - U.^3);

end
