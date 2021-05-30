function [Udot] = randd_slow(t,U)
  % Computes the slow rhs function for the reaction-diffusion problem
  % with fixed linearization
  
  global gamma
  Udot = gamma*(U.^2 - U.^3);
end
