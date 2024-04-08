
% Continuation method using natural continuation
% G(u,nu) = \phi(t,u(nu),nu) - \gamma u(nu) = 0, G: R^n x R -> R^n
% The nonlinear system G is solved using Newton's method
% and the linear system is solved using gmres


% The output file name
filename = ['re_vac_sln_ds_m1',int2str(today)]


% computational parameters for continuation
Nmax = 5;           % maximum Newton iteration
Newt_tol = 1.0e-8;  %  tolerance for Newton iteration
n_eigs=12;          % number of eigenvalues to compute

kmax = 50;  % number of steps along the solution curve
ds = 0.25;  % parameter increments along solution branch


% parameters needed for the time stepper of the sheared electroconvection
K  = 32;                % highest Fourier wave number
Nc  = 24;               % highest power of the Chebyshev poly
%  other parameters defined in SimCode files 
%     (assigned values here only for recording purposes; 
%        changes in these parameters must be made in SimCode files)
rr = 0.56;               % aspect ratio
Re = 0.231; 		% rotational Reynolds number
Pr = 75.8;              % Prandlt number

dt = 2.0e-4;            % time step


% degrees of freedom:    n = 2*(Nc+1)*(2*K+1)
% X_init <--> [u0;omega_init;period_init] :  n+2 x 1
%  Ra (scalar) is continuation parameter

%  Get initial solution (or guess) -> X_init
k_st=1;  % column in which known solution is stored
load ('mat_files/vac_cont_data1.mat','X')  % already computed solutions

% if we have more than one known solution use a secant to get the first guess
%    ds_last = X(end,k_st)-X(end,k_st-1)
%    X_init = X(1:end-1,k_st) + ( X(1:end-1,k_st)-X(1:end-1,k_st-1) ) * ds/ds_last;
%    Ra = X(end,k_st) + ds  

% otherwise
    X_init=X(1:end-1,k_st);
    Ra=X(end,k_st) % + ds

%  possibly make a better initial guess by integrating
%    U_t = X_init(1:end-2);
%    perd = X_init(end)
%     [U_t,Uout] = TS_3x_trunc([0 3*perd],U_t,Ra,dt,K,Nc);
%
%    X_init(1:end-2)=U_t;
%
%   'done'
%


%  initialization of results and diagnostic variables
N = length(X_init); 
X = zeros(N+1,kmax);  % X <--> [u;omega;period;Ra]
mu = []; %zeros(n_eigs,kmax);
V = zeros(N-2,n_eigs,kmax);
eig_conv = zeros(n_eigs,kmax);
time_count_eig = zeros(n_eigs,kmax);
period      = zeros(1,kmax);
phase_speed = zeros(1,kmax);
NEW_res     = zeros(Nmax,kmax);
NEW_abs     = zeros(Nmax,kmax);
GMRES_rel   = zeros(Nmax,kmax);
GMRES_flag  = zeros(Nmax,kmax);

DPhi=zeros(N);
size(DPhi)

run_info = ['the following parameters values have been used \n' ...
             'tolerance for Newton iterations:  Newt_tol: 1e-8 \n'...
             'GMRES residual tolerance 1e-6 and epsilon of the function evaluation is 1e-4\n'];

% initial guess for the rotation wave speed for saving purposes only          
omega_init = X_init(end-1);
period_init = X_init(end);


for k = 1:kmax
    [root,lambda,vec,DPhi,New_iter,flag_eigs,eig_time,New_res_vec,New_abs_vec,Gmres_rel_vec,Gmres_flag_vec]...
                      = AmpMod_NewtKry(X_init,Ra,n_eigs,Nmax,Newt_tol,dt,DPhi,K,Nc);

    X(:,k)   = [root;Ra];               % solution vector [u;omega;period;Ra]
    mu(:,k) = lambda; % = [mu lambda];  % eigenvalues of the map 
    V(:,:,k) = vec;                     % eigenfunctions
    
    eig_conv(k) = flag_eigs;            % eigs flag checking for convergence of eigenvalues 
    time_count_eig(k)   = eig_time;     % CPU time for each computation 
    NEW_res(:,k)      = New_res_vec;    % Newton residual 
    NEW_abs(:,k)      = New_abs_vec;    % Newton absolute error
    GMRES_rel(:,k)    = Gmres_rel_vec;  % GMRES relative error
    GMRES_flag(:,k)   = Gmres_flag_vec; % GMRES flag
    
    period(k) = root(end);
    phase_speed(k) = root(end-1);

    % saving data in a file
    save(filename,'X','mu','V','period','phase_speed','eig_conv','time_count_eig','NEW*','GMRES*',...
        'run_info','omega_init','period_init','ds','kmax','Pr','Re','rr','Nc','K')


    % Initializng the guess for the next point on the solution branch
    Ra = Ra + ds;

    if (k<2)
%  if k=1 and solution is stable can use integration to make a better guess
%      u0 = root(1:end-2);
%      omega = root(end-1)
%      t_end = root(end)

%      'making new initial guess'  
%      [unew, p_full] =  TS_3x_trunc(2*[0 t_end],u0,Ra,dt,K,Nc);
%      'done'
%
%       X_init = [unew;omega;t_end];

% otherwise (for k=1) just use the previous (first) point on the solution curve

      X_init = root;

% if it's not the first point on the curve, use a secant to get the new guess    
    else
      X_init = X(1:end-1,k) + ( X(1:end-1,k)-X(1:end-1,k-1) ) ;

    end

    
end

