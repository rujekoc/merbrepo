% USAGE : driver_randd_merk_fixed
% Run with fixed linearization

% Set up reaction diffusion problem
addpath('../testdescriptions')
addpath('../utility')
[Y0,t0,tf] = randd_setup();
fn  = @randd;
Jn  = @randd_Jacobian;
A   = randd_JacobianF(t0,Y0);
fs  = @randd_slow;

% time parameters
n       = 11;
tout    = linspace(t0,tf,n);                       % intermediate times for solution
hmax    = tout(2)-tout(1);                         % largest time step
h       = hmax*0.5.^(0:6);                         % large time steps

% Initialize problem variables/ allocate space
Y          = zeros(length(Y0),n);
Y(:,1)     = Y0;
err_max    = zeros(1,length(h));
err_rms    = zeros(1,length(h));
max_rate   = zeros(1,length(h)-1);
rms_rate   = zeros(1,length(h)-1);
nffast     = zeros(1,length(h));
nfslow     = zeros(1,length(h));
time       = zeros(1,length(h));

% Get reference solution
tic
hmin = 1e-6;
hmax = 1.0;
opts = odeset('RelTol',1e-13, 'AbsTol',1e-14*ones(size(Y0)),'InitialStep',hmin/10, 'MaxStep',hmax, 'Jacobian',Jn);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
q       = matfile('randd_reference','Writable',true);
Ytrue   = Ytrue';
q.Ytrue = Ytrue;
q.tout  = t;
toc

% % Load reference solution if previously stored
% q = matfile('randd_reference.mat');
% Ytrue = q.Ytrue;

% Solver info
solvers = {@solve_MERK3, @solve_MERK4, @solve_MERK5};
stagesolvers = {'ERK-3-3', 'ERK-4-4', 'ARK5(4)8L[2]SA-ERK'};
stepsolvers = stagesolvers;
snames  = {'MERK3', 'MERK4','MERK5'};
mvals = [10,10,5];
orders  = [3,4,5];

% march over different solvers
for isol = 1:length(solvers)

  % Reset function calls
  nffast     = zeros(1,length(h));
  nfslow     = zeros(1,length(h));

  % Solver information
  solver  = solvers{isol};
  sname   = snames{isol};
  order   = orders(isol);
  rkstage = stagesolvers{isol};
  rkstep  = stepsolvers{isol};
  m       = mvals(isol);

  fprintf('\nRunning convergence test with %s integrator (theoretical order %i)\n',...
  sname,order)
  B = butcher(rkstage);  s = numel(B(1,:))-1;
  D = butcher(rkstep);  s1 = numel(D(1,:))-1;
  fprintf('\nRunning with inner ERK integrators: %s (order = %i)  and %s (order = %i)\n',...
  rkstage,B(s+1,1),rkstep,D(s1+1,1));
  fprintf('\nRunning with m = %i \n',m);

  % March over h values
  for j = 1:length(h)

    % Set timesteps
    hs = h(j);
    hf = h(j)/m;

    % March over tout values
    tstart = tic;
    for k = 2:n
      % Compute solution using MERK solver
      Y0step     = Y(:,k-1);
      [Yout,fastcalls,slowcalls] = solver(A,fs,rkstage,rkstep,Y0step,m,[tout(k-1),tout(k)],hs);
      Y(:,k) = Yout;
      % Update function calls
      nffast(j) = nffast(j)+ fastcalls;
      nfslow(j) = nfslow(j)+ slowcalls;
    end
    telapsed = toc(tstart);
    time(j) = telapsed;

    % Function calls
    fprintf('nffast = %g ,  nfslow = %g .\n',nffast(j),nfslow(j));

    % Error calculation
    err_max(j) = max(max(abs(Y - Ytrue)));
    err_rms(j) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
    fprintf('Accuracy/Work Results, h = %g  , m = %i:\n',hs,m)
    fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max(j), err_rms(j));

    % Rate of convergence calculation
    if j > 1
      max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
      rms_rate(j-1) = log(err_rms(j-1)/err_rms(j))/log(h(j-1)/h(j));
      fprintf('  maxrate = %.5e, rmsrate = %.5e \n',max_rate(j-1),rms_rate(j-1));
    end
    fprintf('Time elapsed = %g.\n',telapsed);

  end % end for hslow values

  % Best-fit convergence rate
  jmax = length(h);
  for j = 2:length(h)
    if (err_max(j)>err_max(j-1)/2)
      jmax = j-1;
      break;
    end
  end
  p1 = polyfit(log(h(1:jmax)),log(err_max(1:jmax)),1);
  fprintf('best-fit convergence max rate = %g\n', p1(1));

  % Best-fit convergence rate
  jmax = length(h);
  for j = 2:length(h)
    if (err_rms(j)>err_rms(j-1)/2)
      jmax = j-1;
      break;
    end
  end
  p2 = polyfit(log(h(1:jmax)),log(err_rms(1:jmax)),1);
  fprintf('best-fit convergence rms rate = %g\n', p2(1));

end % end for solvers
