function [Dni] = randd_Dni(t,dt,U,V)
  % Computes Dni functions for MERB methods for the reaction-diffusion test problem
  %   Dni = F(t+dt,V) - J*V - Nn
  % where
  %   J   = F_U(t,U)
  %   Nn  = F(t,U) - J*U

  global gamma
  global A

  Dni = gamma*( (V - V.^2 - 2*U + 3*U.^2).*V + (U - 2*U.^2).*U );

end
