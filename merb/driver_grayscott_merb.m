function driver_grayscott_merb(isol)
  % driver_grayscott_merb
  % driver to test MERB methods with implicit inner integrators on Gray-Scott
  % problem
  % Input : index for MERB solver
  %       1(no input) - solve_MERB3
  %       2           - solve_MERB4
  %       3           - solve_MERB5
  %       4           - solve_MERB6
  %
  % Rujeko Chinomona
  % Temple University
  % May 2022

  addpath('../testdescriptions')
  addpath('../utility')

  % if no input set isol to 1
  if nargin < 1, isol=1; end

  % Set up gray scott problem
  [Y0,t0,tf] = gray_scott_setup();
  fn = @gray_scott;  % full rhs
  Jfy = @gray_scott_Jacobian; % Jacobian of rhs
  [nrow,ncol] = size(Y0); % size of solution vector
  Jft = @(t,y) sparse(nrow,ncol); % Time derivative
  Nn = @gray_scott_nonlinearity; % rhs nonlinear function
  Dn = @gray_scott_Dni; % useful rhs function

  % Time parameters
  n = 11; % solution outputs
  tout = linspace(t0,tf,n); % intermediate times for output
  hmax = 0.01; % largest slow step size
  if (hmax-(tout(2)-tout(1))>0)
    fprintf('Warning: hmax is bigger than spacing between output times. Pick smaller hmax\n');
    hmax = tout(2) - tout(1);
  end
  h = hmax*0.5.^(0:7); % slow step sizes for convergence tests
  mvals = 10; % ratio between slow and fast step sizes

  % Get reference solution
  fprintf('Compute reference solution with ode15s ... \n')
  hmin = 1e-6;
  opts = odeset('RelTol',1e-13, 'AbsTol',1e-14*ones(size(Y0)),'InitialStep',hmin/10, 'MaxStep',hmax, 'Jacobian',Jfy);
  [t,Ytrue] = ode15s(fn, tout, Y0, opts);

  % Save results to .mat files
  q       = matfile('gscott_reference','Writable',true);
  q.dateandtime = now;
  q.Ytrue = Ytrue;
  q.tout = t;

  % % Load reference solution if previously stored
  % q = matfile('gscott_reference.mat');
  % Ytrue = q.Ytrue;

  % Initialize problem variables/allocate space
  Y          = zeros(length(Y0),n);
  Y(:,1)     = Y0;
  err_max    = zeros(length(mvals),length(h));
  err_rms    = zeros(length(mvals),length(h));
  max_rate   = zeros(length(mvals),length(h)-1);
  rms_rate   = zeros(length(mvals),length(h)-1);
  nffast     = zeros(length(mvals),length(h));
  nfslow     = zeros(length(mvals),length(h));
  time       = zeros(length(mvals),length(h));

  % Solver info
  solvers = {@solve_MERB3, @solve_MERB4, @solve_MERB5,@solve_MERB6};
  stepsolvers = {'RadauIIA-2-3-IRK','LobattoIIIC-3-4-IRK', ...
  'RadauIIA-3-5-IRK','LobattoIIIC-4-6-IRK'};
  stagesolvers = stepsolvers;   % inner stages solved at same accuracy as step solution
  snames  = {'MERB3','MERB4','MERB5','MERB6'}; % method names
  orders  = [3,4,5,6]; % Order of convergence

  % Set Solver information
  solver  = solvers{isol};
  sname   = snames{isol};
  order   = orders(isol);
  rkstage = stagesolvers{isol};
  rkstep  = stepsolvers{isol};

  fprintf('\nRunning convergence test with %s integrator (theoretical order %i)\n',...
  sname,order)
  Btable = butcher(rkstage);  s = numel(Btable(1,:))-1;
  Dtable = butcher(rkstep);  s1 = numel(Dtable(1,:))-1;
  fprintf('\nRunning with inner RK integrators: %s (order = %i)  and %s (order = %i)\n',...
  rkstage,Btable(s+1,1),rkstep,Dtable(s1+1,1));

  for i=1:length(mvals)
    m = mvals(i);
    fprintf('\nRunning with m = %i \n',m);

    % March over h values
    for j = 1:length(h)

      % Set timesteps
      hs = h(j);

      % March over tout values
      tstart = tic;
      for k = 2:n
        % Compute solution using MERB solver
        Y0step     = Y(:,k-1);
        [Yout,fastcalls,slowcalls]= solver(fn,Jfy,Jft,rkstage,rkstep,Y0step,m,[tout(k-1),tout(k)],hs,true,Nn,Dn);
        Y(:,k) = Yout;
        % Update function calls
        nffast(i,j) = nffast(i,j)+ fastcalls;
        nfslow(i,j) = nfslow(i,j)+ slowcalls;
      end
      telapsed = toc(tstart);
      time(i,j) = telapsed;

      % Error calculation
      err_max(i,j) = max(max(abs(Y' - Ytrue)));
      err_rms(i,j) = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
      fprintf('Accuracy/Work Results, h = %g  , m = %i:\n',hs,m)
      fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max(i,j), err_rms(i,j));

      % Rate of convergence calculation
      if j > 1
        max_rate(i,j-1) = log(err_max(i,j-1)/err_max(i,j))/log(h(j-1)/h(j));
        rms_rate(i,j-1) = log(err_rms(i,j-1)/err_rms(i,j))/log(h(j-1)/h(j));
        fprintf('  maxrate = %.5e, rmsrate = %.5e \n',max_rate(i,j-1),rms_rate(i,j-1));
      end

      % Function calls
      fprintf('nffast = %g ,  nfslow = %g .\n',nffast(i,j),nfslow(i,j));
      fprintf('Time elapsed = %g.\n',telapsed);

    end % end for hslow values

    % Best-fit convergence rate
    jmax = length(h);
    for j = 2:length(h)
      if (err_max(i,j)>err_max(i,j-1))
        jmax = j-1;
        break;
      end
    end
    p1 = polyfit(log(h(1:jmax)),log(err_max(i,1:jmax)),1);
    fprintf('best-fit convergence max rate = %g\n', p1(1));

    % Best-fit convergence rate
    jmax = length(h);
    for j = 2:length(h)
      if (err_rms(i,j)>err_rms(i,j-1)/2)
        jmax = j-1;
        break;
      end
    end
    p2 = polyfit(log(h(1:jmax)),log(err_rms(i,1:jmax)),1);
    fprintf('best-fit convergence rms rate = %g\n', p2(1));

  end
end
