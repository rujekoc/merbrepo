function [Unew,nffast,nfslow] = solve_MERB6(Fn,Jfy,Jft,stagesolver,stepsolver, U0, m, tinterval,H,implicitsolve,N,Dn)
  % usage: [Unew,nffast,nfslow] = solve_MERB6(Fn,Jfy,Jft,stagesolver,stepsolver, U0, m, tinterval,H,implicitsolve,N,Dn)
  %
  % Sixth order multirate exponential Rosenbrock solver for the vector-valued ODE problem
  %     u' = F(t,u), t in tinterval, y in R^m,
  %     u(t0) = [u1(t0), u2(t0), ..., um(t0)]'.
  %
  % Inputs:
  %     Fn          = function handle for F(t,u(t))
  %     Jfy         = function handle for (Jacobian) = F_u(t,u(t))
  %     Jft         = function handle for F_t(t,u(t))
  %     stagesolver = string holding function name for  RK method
  %                 Butcher table formatting
  %                 B = [ci A;
  %                      q b]
  %                 Here, ci is a vector of stage time fractions (s-by-1),
  %                    A is a matrix of Butcher coefficients (s-by-s),
  %                    q is an integer denoting the method order of accuracy,
  %                    b is a vector of solution weights (1-by-s),
  %     stepsolver  = string holding function name for RK method
  %     tinterval   = [t_n,t_n+1] initial and final time
  %     U0          = initial value array (column vector of length m)
  %     H           = slow (macro) time step
  %   implicitsolve = true/logical 1 if using implicit RK for inner solve,
  %                   false/ logical 0 for explicit RK
  % Optional inputs:
  %     N           = function handle for nonlinearity in F: F-F_u*u
  %     Dn          = function handle for F(t+dt,V) - F_U(t,U)*V - F(t,u)-F_U(t,U)*U - dt*F_t(t,U)
  %
  % Outputs:
  %     Unew    = column vector with solution at final time
  %     nffast = number of fast function calls
  %     nfslow = number of slow function calls
  %
  %
  % Rujeko Chinomona
  % Department of Mathematics
  % Temple University
  % February 2022

  % flags for definition of N and Dn functions
  user_N  = false;
  user_Dn = false;
  if (nargin > 10)
    user_N = true;
  end
  if (nargin > 11)
    user_Dn = true;
  end

  % Info from butcher table for inner ODE solvers
  B = butcher(stagesolver);
  D = butcher(stepsolver);

  % Set problem parameters
  c       = [0,1/9,1/10,1/7,1/10,1/9,1/8];      % ExpRB stage abscissae
  n       = 0;
  Unew    = U0;                                 % initial condition for this step
  tn      = tinterval(1);                       % initial time
  ONEMSM  = 1-sqrt(eps);                        % coefficient to account for floating-point roundoff
  nffast  = 0;                                  % initialize function evaluations
  nfslow  = 0;
  rtol    = 1e20;
  atol    = rtol*ones(size(U0));

  while tn < tinterval(2)*ONEMSM

    H  = min([H,tinterval(2)-tn]);

    % Set Un1
    Y0  = Unew;

    % Determine micro time step
    hf = c(3)*H/ceil(c(3)*m);
    tsteps = [0,c(3)*H,c(2)*H];

    % Set Jacobian information and slow function
    Jac     = Jfy(tn,Unew);
    JacT    = Jft(tn,Unew);
    if (user_N)
      Nn    = N(tn,Unew);
    else
      Nn    = Fn(tn,Unew) - Jac*Unew;
    end
    nfslow  = nfslow + 1;               % Update number of slow function calls

    % Define slow contribution function
    if (user_Dn)
      Dni   = @(t,c,y) Dn(t,c*H,Unew,y);
    else
      Dni   = @(t,c,y) Fn(t+c*H,y) - Jac*y - Nn  - c*H*JacT;
    end

    % Set up right hand side for modified ODE to solve for 2nd stage
    pn2 = @(t)   Nn + t*JacT;
    fcn = @(t,y) Jac*y + pn2(t);
    Jfcn = @(t,y) Jac;

    %Solve for 2nd and 3rd stage
    if implicitsolve
      [~,Y,nsteps,lits,~] = solve_IRK(fcn,Jfcn,tsteps,Y0,B,rtol,atol,hf,hf,hf);
      nflocal = (nsteps + lits)*(numel(B(1,:))-1);
    else
      [~,Y,nflocal] = solve_ERKfast(fcn,tsteps,Y0,B,hf);
    end
    Z2 = Y(:,end);
    Z3 = Y(:,2);
    nffast = nffast + nflocal;          % Update number of fast function calls

    % Solve for slow function contribution
    Dn2 = Dni(tn,c(2),Z2);
    Dn3 = Dni(tn,c(3),Z3);
    nfslow   = nfslow + 2;              % Update number of slow function calls

    % Determine micro time step
    hf = c(5)*H/ceil(c(5)*m);
    tsteps = [0,c(5)*H,c(6)*H,c(7)*H,c(4)*H];

    % Set up right hand side for modified ODE to solve for 3rd and 4th stages
    pn4 = @(t) Nn + ...
          (t/H)*(H*JacT + ...
                 (t/H)*((c(3)*Dn2/c(2)^2/(c(3) - c(2)) + c(2)*Dn3/c(3)^2/(c(2)-c(3))) + ...
                        (t/H)*(-Dn2/c(2)^2/(c(3) - c(2)) - Dn3/c(3)^2/(c(2) - c(3)))));
    fcn  = @(t,y) Jac*y + pn4(t);

    %Solve for Un3 and Un4
    if implicitsolve
      [~,Y,nsteps,lits,~] = solve_IRK(fcn,Jfcn,tsteps,Y0,B,rtol,atol,hf,hf,hf);
      nflocal = (nsteps + lits)*(numel(B(1,:))-1);
    else
      [~,Y,nflocal] = solve_ERKfast(fcn,tsteps,Y0,B,hf);
    end
    Z5 = Y(:,2);
    Z6 = Y(:,3);
    Z7 = Y(:,4);
    Z4 = Y(:,end);
    nffast = nffast + nflocal;          % Update number of fast function calls

    % Solve for slow function contribution
    Dn4 = Dni(tn,c(4),Z4);
    Dn5 = Dni(tn,c(5),Z5);
    Dn6 = Dni(tn,c(6),Z6);
    Dn7 = Dni(tn,c(7),Z7);
    nfslow   = nfslow + 4;                 % Update number of slow function calls

    g = zeros(1,7);
    a = zeros(1,7);
    b = zeros(1,7);
    e = zeros(1,7);

    for i = 4:7
      g(i) = 1/c(i)^2;
      for k = 4:7
        if k ~= i
          g(i) = g(i)/(c(i) - c(k));
        end
      end
    end

    for i = 4:7
      a(i) = g(i);
      for k = 4:7
        if k~=i
          a(i) = a(i)*c(k);
          b(i) = b(i) + c(k);
        end
      end
      b(i) = b(i)*g(i);
    end

    e(4) = (c(5)*c(6) + c(6)*c(7) + c(5)*c(7))*g(4);
    e(5) = (c(4)*c(6) + c(6)*c(7) + c(4)*c(7))*g(5);
    e(6) = (c(4)*c(5) + c(5)*c(7) + c(4)*c(7))*g(6);
    e(7) = (c(4)*c(5) + c(5)*c(6) + c(4)*c(6))*g(7);

    % Set up right hand side for modified ODE to solve for u_{n+1}
    qn = @(t) Nn + ...
         (t/H)*(H*JacT + ...
                (t/H)*(-a(4)*Dn4 - a(5)*Dn5 - a(6)*Dn6 - a(7)*Dn7 + ...
                       (t/H)*(e(4)*Dn4 + e(5)*Dn5 + e(6)*Dn6 + e(7)*Dn7 + ...
                              (t/H)*(-b(4)*Dn4 - b(5)*Dn5 - b(6)*Dn6 - b(7)*Dn7 + ...
                                     (t/H)*(g(4)*Dn4 + g(5)*Dn5 + g(6)*Dn6 + g(7)*Dn7)))));
    fcn = @(t,y) Jac*y + qn(t);

    % Set up micro time step
    hf = H/m;

    % Solve for u_{n+1} on [0,H]
    if implicitsolve
      [~,Y,nsteps,lits,~] = solve_IRK(fcn,Jfcn,[0,H],Y0,D,rtol,atol,hf,hf,hf);
      nflocal = (nsteps+lits)*(numel(D(1,:))-1);
    else
      [~,Y,nflocal] = solve_ERKfast(fcn,[0,H],Y0,D,hf);
    end
    Unew = Y(:,end);
    nffast = nffast + nflocal;            % Update number of fast function calls

    % Update time step
    tn = tn + H;
    n = n+1;

  end  %end while
end    %end function
