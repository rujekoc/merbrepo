% USAGE : driver_randd_mrigark_fixed
% Does fixed linearization

% Set up reaction diffusion problem
addpath('../testdescriptions')
addpath('../utility')

[Y0,t0,tf] = randd_setup();
fn  = @randd;
ff  = @randd_fast;
fs  = @randd_slow;
Jn  = @randd_Jacobian;

% time parameters
n       = 11;
tout    = linspace(t0,tf,n);                       % intermediate times for solution
hmax    = tout(2)-tout(1);                         % largest time step
h       = hmax*0.5.^(0:6);                         % large time steps

% Get reference solution
tic
hmin = 1e-6;
hmax = 1.0;
opts = odeset('RelTol',1e-13, 'AbsTol',1e-14*ones(size(Y0)),'InitialStep',hmin/10, 'MaxStep',hmax, 'Jacobian',Jn);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
Ytrue     = Ytrue';
q         = matfile('randd_reference','Writable',true);
q.Ytrue   = Ytrue;
q.tout    = t;
toc

% Solver approach 1: exp/exp (4th+4th order)
disp('Solver approach 1: explicit+explicit (4th+4th order)')
BoName = 'MRI-GARK-ERK45a';
co = [0; 1/5; 2/5; 3/5; 4/5; 1];
Gammas{1} = [1/5 0 0 0 0 0; ...
             -53/16, 281/80, 0 0 0 0; ...
             -36562993/71394880, 34903117/17848720, -88770499/71394880, 0 0 0; ...
             -7631593/71394880, -166232021/35697440, 6068517/1519040, 8644289/8924360, 0 0; ...
             277061/303808, -209323/1139280, -1360217/1139280, -148789/56964, 147889/45120 0];
Gammas{2} = [0 0 0 0 0 0; ...
             503/80, -503/80, 0 0 0 0; ...
             -1365537/35697440, 4963773/7139488, -1465833/2231090, 0 0 0; ...
             66974357/35697440, 21445367/7139488, -3, -8388609/4462180, 0 0; ...
             -18227/7520, 2, 1, 5, -41933/7520 0];
BiName = 'ERK-4-4';
Bi = butcher(BiName);
mvals = 1;
[err_max,err_rms,max_rate,rms_rate,time,nffast,nfslow] = ...
dotest_fixed(co,Gammas,Bi,fs,ff,tout,Y0,h,mvals,Ytrue);

% Clear out tables for next run
clear co Gammas Bi mvals

% Solver approach 2: exp/exp (3rd+3rd order)
disp('Solver approach 2: explicit+explicit (3rd+3rd order)')
BoName = 'MRI-GARK-ERK33a';
co = [0; 1/3; 2/3; 1];
Gammas{1} = [1/3 , 0 , 0 , 0; ...
            -1/3 , 2/3 , 0 , 0; ...
              0 , -2/3 , 1, 0];
Gammas{2} = [ 0 , 0 , 0 , 0; ...
              0 , 0 , 0 , 0; ...
              1/2 , 0 , -1/2, 0];
BiName = 'ERK-3-3';
Bi = butcher(BiName);
mvals = 5;
[err_max,err_rms,max_rate,rms_rate,time,nffast,nfslow] = ...
dotest_fixed(co,Gammas,Bi,fs,ff,tout,Y0,h,mvals,Ytrue);

% Clear out tables for next run
clear co Gammas Bi mvals
