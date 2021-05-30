function [Udot] = randd_fast(t,U)
  % Computes the fast linear rhs function for the reaction-diffusion test problem
  % with fixed linearization

  global A
  Udot = A*U;
end
