function [J] = randd_Jacobian(t,U)
  % Computes the full linearization of the rhs of the reaction-diffusion test problem
  % with dynamic linearization
  
  global gamma
  global A

  J = A + gamma*diag(2*U - 3*U.^2);

end
