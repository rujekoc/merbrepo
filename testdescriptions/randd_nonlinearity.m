function [Nn] = randd_nonlinearity(t,U)
  % Computes the nonlinear portion left after full linearization of rhs
  %   Nn = F(t,U) - J(t,U)*U

  global gamma
  Nn = gamma*( -U.^2 + 2*U.^3 );

end
