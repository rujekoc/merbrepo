function [Nn] = randd_nonlinearity_gen(tn,Un)
  % Generates the Nn function for the reaction-diffusion test problem
  % in the case of dynamic linearization for MERK and MRIGARK methods.
  %   Nn = F(t,U) - J(t,U)*U

  global gamma
  Nn = @(t,U) gamma*(U.^2 - U.^3 - 2*Un.*U + 3*Un.^2.*U);

end
