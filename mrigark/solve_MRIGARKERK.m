function [tvals,Y,ns,nf,ierr] = solve_MRIGARKERK(fs,ff,tvals,Y0,co,G,Bi,hs,hf,Nn,Jfy)
  % usage: [tvals,Y,ns,nf,ierr] = solve_MRIGARKERK(fs,ff,tvals,Y0,co,G,Bi,hs,hf,Nn,Jfy)
  % Can be used for both dynamic or fixed linearization approaches
  % Fixed time step Multirate Infinitesimal Step (MRI) solver for the vector-valued ODE problem
  %     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
  %     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
  % This solver is designed to handle only EXPLICIT fast and slow methods
  % The individual time steps are performed using the step_MRI function;
  % that in turn iterates over the stages for the slow-time-scale method,
  % and calls structure-specific routines for each stage.
  %
  % Inputs:
  %     fs     = function handle for (slow) ODE RHS (fixed linearization)
  %     ff     = function handle for (fast) ODE RHS (fixed linearization)
  %     tvals  = array of desired output times, [t0, t1, t2, ..., tN]
  %     Y0     = solution vector at start of step (column vector of length n)
  %     co     = Outer 'slow' method abcissae. Here we assume the following
  %                    0 < co(1) < co(2) < ... < co(end) <= 1
  %               so there is no division by zero
  %     G      = cell array of MRI "Gamma" matrices, {G0, G1, ...}
  %     Bi     = Butcher table for a single step of an 'inner' (fast) method
  %              (can be ERK, DIRK or IRK)
  %                 Bi = [ci Ai;
  %                       qi bi;
  %                       pi di ]
  %              All components have the same role as with Bi; we
  %              assume that Bi encodes a method with si stages
  %     hs     = step size to use for slow time scale
  %     hf     = desired step size to use for fast time scale,
  %                  hf <= hs*min_{i}(co(i+1)-co(i))
  %              Note: this is only a target step size; in fact we
  %              will determine each substepping interval and find
  %              hinner <= hi such that we take an integer number of
  %              steps to subcycle up to each outer stage time.
  %    Optional inputs for dynamic linearization:
  %     Nn     = function handle for nonlinearity generating function
  %     Jfy    = function handle for Jacobian of fn
  % Outputs:
  %     tvals  = the same as the input array tvals
  %     Y      = [y1(t0+hs), y2(t0+hs), ..., yn(t0+hs)].
  %     ns     = number of 'slow' time steps taken by method
  %     nf     = number of 'fast' function calls taken by method
  %     ierr   = flag denoting success (0) or failure (1)
  %
  % Daniel R. Reynolds, Rujeko Chinomona
  % Department of Mathematics
  % Southern Methodist University
  % May 2021
  % All Rights Reserved

  % flags for dynamic linearization
  dynamic  = false;
  if (nargin > 9)
    dynamic = true;
  end
  if (nargin == 10)
    error('Total number of inputs is 9 or 11')
    return;
  end

  % initialize output arrays
  N = length(tvals)-1;
  n = length(Y0);
  Y = zeros(n,N+1);
  Y(:,1) = Y0;

  % initialize diagnostics and error flag
  ns = 0;
  nf = 0;
  ierr = 0;

  % set the solver parameters
  ONEMSM = 1-sqrt(eps);  % coefficient to account for floating-point roundoff

  % initialize temporary variables
  t = tvals(1);
  Ynew = Y0;

  % iterate over output time steps
  for tstep = 1:N

    % loop over internal time steps to get to desired output time
    while (t < tvals(tstep+1)*ONEMSM)

      % bound internal time step
      h = min([hs, tvals(tstep+1)-t]);   % stop at output times

      % Check if implementing dynamic linearization
      if (dynamic)
        % Set Jacobian
        A = Jfy(t,Ynew);
        % Set nonlinear function
        fs = Nn(t,Ynew);
        % Set fast function
        ff = @(tau,y) A*y;
      end

      % call MIS stepper to do the work, increment counters
      [Ynew,m,ierr] = step_MRI(fs,ff,t,Ynew,co,G,Bi,h,hf);
      ns    = ns + 1;
      nf    = nf + m;

      % handle error flag
      if (ierr ~= 0)
        fprintf('Error: solve_MRI encountered an error at t = %g, aborting.\n',t);
        return;
      end

      % update current time
      t = t + h;

    end

    % store updated solution in output array
    Y(:,tstep+1) = Ynew;

  end  % time output loop

  % end function
end


function [Y,m,ierr] = step_MRI(fs,ff,t0,Y0,co,G,Bi,hs,hf)
  % usage: [Y,m,ierr] = step_MRI(fs,ff,t0,Y0,co,G,Bi,hs,hf)
  %
  % This routine performs a single step of an explicit-at-slow multirate
  % infinitesimal step (MIS) method for the vector-valued ODE problem
  %     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
  %     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
  %
  % Inputs:
  %     fs     = function handle for (slow) ODE RHS
  %     ff     = function handle for (fast) ODE RHS
  %     t0     = value of independent variable at start of step
  %     Y0     = solution vector at start of step (column vector of length n)
  %     co     = Outer 'slow' method abcissae.  The MRI method requires that
  %              these be sorted, i.e.
  %                    0 <= co(1) <= co(2) <= ... <= co(end) <= 1
  %              and (potentially) padded such that co(end)=1
  %     G      = cell array of MRI "Gamma" matrices, {G0, G1, ...}
  %     Bi     = Butcher table for a single step of an 'inner' (fast) method
  %              (can be ERK, DIRK or IRK)
  %                 Bi = [ci Ai;
  %                       qi bi;
  %                       pi di ]
  %              All components have the same role as with Bi; we
  %              assume that Bi encodes a method with si stages
  %     hs     = step size to use for slow time scale
  %     hf     = desired step size to use for fast time scale,
  %                  hf <= hs*min_{i}(co(i+1)-co(i))
  %              Note: this is only a target step size; in fact we
  %              will determine each substepping interval and find
  %              hinner <= hi such that we take an integer number of
  %              steps to subcycle up to each outer stage time.
  %
  % Outputs:
  %     Y     = [y1(t0+hs), y2(t0+hs), ..., yn(t0+hs)].
  %     m     = actual number of fast function calls.
  %     ierr  = flag denoting success (0) or failure (1)

  % check that abscissae satisfy MRI assumptions
  so = length(co);          % number of stages
  Delta_co = co(2:end)-co(1:end-1);
  nGammas = length(G);
  if ( (abs(co(1)) > 100*eps) || (min(Delta_co) < -100*eps) || (abs(co(so)-1) > 100*eps) )
    error('Error: co does not satisfy MRI assumptions')
  end

  % initialize outputs
  m = 0;
  n = length(Y0);
  Y = reshape(Y0,n,1);
  z = Y;
  tcur = t0;
  t = t0;
  ierr = 0;

  % initialize temporary variables
  Fs = zeros(n,length(co));        % storage for all slow RHS values

  % first outer stage is always explicit
  Fs(:,1) = fs(t0,Y);

  %   time interval
  tspan = [0,hs];

  % iterate over remaining outer stages
  for stage = 2:length(co)
    % number of  internal time steps
    ni = ceil(Delta_co(stage-1)*hs/hf);

    %   step size
    hi = hs/ni;

    % utility routine to construct column vector with powers of (theta/H)
    thetapow = @(theta) ((theta/hs).^(0:nGammas-1))';

    % Gamma matrix for this fast RHS
    Gamma = zeros(so,nGammas);
    for col=1:nGammas
       for row=1:stage
          Gamma(row,col) = G{col}(stage-1,row);
       end
    end

    fi = @(theta,v) Delta_co(stage-1)*ff(tcur + Delta_co(stage-1)*theta,v) + Fs*Gamma*thetapow(theta);

    % call inner RK method solver to perform substepping
    % ERK inner method
    [~,V,mi] = solve_ERKfast(fi,tspan,z,Bi,hi);
    m = m + mi;

    z = V(:,end);

    % construct new slow stage RHS and return
    tcur = t + co(stage)*hs;
    Fs(:,stage) = fs(tcur,z);
  end
  Y = z;

end
