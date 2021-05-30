function [err_max,err_rms,max_rate,rms_rate,time,nffast,nfslow] = dotest_dynamic(co,Gammas,Bi,fs,ff,Nn,Jn,tvals,Y0,h,mvals,Ytrue)
  % USAGE:[err_max,err_rms,max_rate,rms_rate,time,nffast,nfslow] = dotest_dynamic(co,Gammas,Bi,fs,ff,Nn,Jn,tvals,Y0,h,mvals,Ytrue)
  % Calls solve_MRIGARKERK in the dynamic linearization setting and outputs
  % results
  %
  % INPUT: co, Gammas - coefficients of the MRI-GARK method
  %                Bi - coefficients for explicit RK method to compute modified ODEs
  %                fs - function handle for slow ODE rhs
  %                ff - function handle for fast ODE rhs
  %                Nn - function handle for nonlinear portion of problem
  %                Jn - Jacobian of ODE rhs
  %             tvals - solution output times
  %               Y0  - initial condition
  %                h  - slow time step sizes
  %            mvals  - integer ratio of slow time step to fast time step
  %             Ytrue - reference solution
  %
  % OUTPUT :
  %           err_max - maximum absolute errors
  %           err_rms - root mean square errors
  %          max_rate - rate of convergence (max errors)
  %          rms_rate - rate of convergence (rms errors)
  %              time - runtime
  %            nffast - fast function evaluations
  %            nfslow - slow function evaluations

  % Initialize problem variables/ allocate space
  Y          = zeros(length(Y0),length(tvals));
  Y(:,1)     = Y0;
  err_max    = zeros(length(mvals),length(h));
  err_rms    = zeros(length(mvals),length(h));
  max_rate   = zeros(length(mvals),length(h)-1);
  rms_rate   = zeros(length(mvals),length(h)-1);
  nffast     = zeros(length(mvals),length(h));
  nfslow     = zeros(length(mvals),length(h));
  time       = zeros(length(mvals),length(h));

  for i=1:length(mvals)
    m = mvals(i);
    fprintf('\nRunning with m = %i \n',m);
    for j=1:length(h)

      % Set timesteps
      hs = h(j);
      hf = hs/m;

      % March over tvals
      tstart = tic;
      [tout,Y,ns,nf,ierr] = solve_MRIGARKERK(fs,ff,tvals,Y0,co,Gammas,Bi,hs,hf,Nn,Jn);
      telapsed = toc(tstart);
      time(i,j) = telapsed;

      if (ierr ~= 0)   % on a solver failure, just skip to the next h
        err_max(i,j) = 1e4;
        err_max(i,j) = 1e4;
        fprintf('    Solver failure, skipping to next h value\n')
        continue
      end

      % Error calculation
      err_max(i,j) = max(max(abs(Y-Ytrue)));
      err_rms(i,j) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
      fprintf('Accuracy/Work Results, h = %g  , m = %i:\n',hs,m)
      fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max(i,j), err_rms(i,j));

      % Rate of convergence calculation
      if (j>1)
        max_rate(i,j-1) = log( err_max(i,j-1)/err_max(i,j) ) / log( h(j-1)/h(j) );
        rms_rate(i,j-1) = log( err_rms(i,j-1)/err_rms(i,j) ) / log( h(j-1)/h(j) );
        fprintf('  maxrate = %.5e, rmsrate = %.5e \n',max_rate(i,j-1),rms_rate(i,j-1));
      end

      % Function calls and steps
      si = numel(Bi(1,:))-1;
      so = length(co)- 1;
      nffast(i,j) = nf;
      nfslow(i,j) = ns*so;
      fprintf('    slow steps = %i, fast steps = %i\n', ns, nf/si);
      fprintf('nffast = %g ,  nfslow = %g .\n',nffast(i,j),nfslow(i,j));

      fprintf('Time elapsed = %g.\n',telapsed);
    end

    % Best-fit convergence rate
    jmax = length(h);
    for j = 2:length(h)
      if (err_max(i,j)>err_max(i,j-1)/2)
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
