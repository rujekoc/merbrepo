function [Unew,nffast,nfslow] = solve_MERB3(Fn,Jfy,Jft,stagesolver,stepsolver, U0, m, tinterval,H,implicitsolve,N,Dn)
  % usage: [Unew,nffast,nfslow] = solve_MERB3(Fn,Jfy,Jft,stagesolver,stepsolver, U0, m, tinterval,H,implicitsolve,N,Dn)
  %
  % Third order multirate exponential Rosenbrock solver for the vector-valued ODE problem
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
  c       = [0,1/2];             % ExpRB stage abscissae
  n       = 0;
  Unew    = U0;                  % initial condition for this step
  tn      = tinterval(1);        % initial time
  ONEMSM  = 1-sqrt(eps);         % coefficient to account for floating-point roundoff
  nffast  = 0;                   % initialize function evaluations
  nfslow  = 0;
  rtol    = 1e20;                % tolerances for implicit RK solver
  atol    = rtol*ones(size(U0));

  while tn < tinterval(2)*ONEMSM

    H  = min([H,tinterval(2)-tn]);

    % Set Un1
    Y0  = Unew;

    % Determine micro time step
    hf = c(2)*H/ceil(c(2)*m);

    % Set Jacobian information and slow function
    Jac     = Jfy(tn,Unew);
    JacT    = Jft(tn,Unew);
    if (user_N)
      Nn    = N(tn,Unew);
    else
      Nn    = Fn(tn,Unew) - Jac*Unew;
    end
    nfslow   = nfslow + 1;             % Update number of slow function calls

    % Define slow contribution function
    if (user_Dn)
      Dni   = @(t,c,y) Dn(t,c*H,Unew,y);
    else
      Dni   = @(t,c,y) Fn(t+c*H,y) - Jac*y - Nn  - c*H*JacT;
    end

    % Set up right hand side for modified ODE to solve for 2nd stage
    pn  = @(t)   Nn + t*JacT;
    fcn = @(t,y) Jac*y + pn(t);
    Jfcn = @(t,y) Jac;

    % Solve for Un2
    if implicitsolve
      [~,Y,nsteps,lits,~] = solve_IRK(fcn,Jfcn,[0,c(2)*H],Y0,B,rtol,atol,hf,hf,hf);
      nflocal = (nsteps + lits)*(numel(B(1,:))-1);
    else
      [~,Y,nflocal] = solve_ERKfast(fcn,[0,c(2)*H],Y0,B,hf);
    end
    Z2 = Y(:,end);
    nffast = nffast + nflocal;         % Update number of fast function calls

    % Solve for slow function contribution
    Dn2 = Dni(tn,c(2),Z2);
    nfslow   = nfslow + 1;             % Update number of slow function calls

    % Set up right hand side for modified ODE to solve for u_{n+1}
    qn  = @(t)   Nn + (t/c(2)/H)^2*Dn2 + t*JacT;
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
    nffast = nffast + nflocal;         % Update number of function calls

    % Update time step
    tn = tn + H;
    n = n+1;

  end %end while
end   %end function
