function [X,mu,V,DMonoD,New_iter,flag_eigs,eig_time,New_res_vec,New_abs_vec,Gmres_rel_vec,Gmres_flag_vec] ...
                        = New_Kry_solver_vac(X,Ra,n_eigs,imax,abs_tol,dt,DMonoD,K,Nc)

% Newton iterations using gmres 

% INPUT:
% X0    mx1 initial point 
% Ra    continuation parameter value
% n_eigs   number of eigenvalues to compute
% abs_tol   tolerance of the error set at 1e-8 by default
% imax  number of iteration set at 5 by default
% dt  time step 
% DPhi -> DMonoD  preconditioning matrix
%  K, Nc  grid size 

% OUTPUT:
% root  1xm root of the function
% mu   1xn_eigs eigenvalues
% V    NxNxn_eigs eigenfunctions
% New_iter  1x1  number of iterations 
% ... and diagnostics
% 

global gmit  % a counter to follow progress of iterations

format shorte

%  default values for computational parameters
if nargin < 5
  abs_tol = 1.0e-8    
    if nargin < 4
      imax = 5
    end
end

imax;

%  initializing variables for computation and diagnostics
New_res_vec = zeros(imax,1);
New_abs_vec = zeros(imax,1);
Gmres_rel_vec = zeros(imax,1);
Gmres_flag_vec = zeros(imax,1);
test = true;
New_iter = 0;
New_res = 1e3;


    u0 = X(1:end-2); % guess of solution (used in auxiliary condition)
    omega = X(end-1) % guess of phase speed
    t_end = X(end)   % guess of period
    tspan = [0 t_end]

    %  compute deriv of u0 wrt theta (for auxiliary condition)
    u_theta_0 = dtheta_electro_trunc(u0,K,Nc);

    %  compute phi = Phi_{period}(u0) for linearization
    unew = u0;
    Nstep=round(t_end/dt);
    dt_now=t_end/Nstep;

    [phi, p_full] =  TS_3x_trunc(tspan,unew,Ra,dt_now,K,Nc);

    % computing initial residual (for reference only)
    gamma = rotation_trunc(unew,omega,tspan,Nc,K);

    b = -[(phi-gamma) ; ...
          u_theta_0.'*(unew); ...
             u0.'*unew - u0.'*u0];

    New_res = norm(b,'inf')


%  load preconditioning matrix, if it exist for current Ra
  in_pmat_fname=['mat_files/PreCond_Rep231_Nc24_K32_Ra',int2str(10*Ra),'_it1_dt_5em4.mat']

if (exist(in_pmat_fname,'file') == 2)
  'loading existing DMonoD'
  load(in_pmat_fname,'DMonoD');
  size(DMonoD)
else
%  'loading first'
%  load('full_matrix_Nc24_K32_Pr18_2.mat','DMonoD')
   'keep same DMonoD'
end


 %%%%%%%%%% start of Newton iteration loop %%%%%%%%%

while (test)

  %%%%  for selecting preconditioner %%%%
  if ( (New_iter == 1) ) & (mod(round(10*Ra),5)==0 )
    
    %%%  computing preconditioner DMonoD %%%

    % dt=5e-4  is set in mat_vec_prod_MonoD  

    'computing preconditioner'
    
    %  output file for preconditioning matrix 
      full_fname=['mat_files/PreCond_Rep231_Nc24_K32_Ra',int2str(10*Ra),'_it1_dt_5em4.mat']

    Mx = @(v) matvec_prod_MonoD(unew,phi,u0,u_theta_0,omega,tspan,Ra,K,Nc,v);


    side_MD=length(b)
    bin_sz = 200
    nbins = 16

    DMonoD=zeros(side_MD);

    %%   Generate DPhi %%

       parfor mmm=bin_sz*nbins+1:side_MD  % replace 'parfor' <-> 'for' if non-parallel
          v=zeros(length(b),1);
          v(mmm)=1;
          DMonoD(:,mmm)=Mx(v);
       end

       save(full_fname,'DMonoD','K','Nc','Ra')


       for kk = 1:nbins
             kk
          parfor mmm=(kk-1)*bin_sz+1:kk*bin_sz; % replace parfor <-> for if non-parallel
             v=zeros(length(b),1);
             v(mmm)=1;
             DMonoD(:,mmm)=Mx(v);
          end

       save(full_fname,'DMonoD','K','Nc','Ra')

       end  
       %%%  end compute DMonoD %%%


    t=2e-4

  else
    Ra
    'keeping same DMonoD as previous'

  end
  %%%  end for loading Preconditioner computation 


  %  define action of linearization (approx with finite diff)

      Ax = @(v) matvec_prod_vac(unew,phi,u0,u_theta_0,omega,tspan,dt,Ra,K,Nc,v);

  %  define preconditioner for gmres function
  M_bs_y = @(v) DMonoD\v;


  gmit=0;
  fprintf('\n gmres iteration:  ')

  %  compute Newton increment (delta_X) using gmres

  [delta_X,flag_gmres,relres,gmres_iter,resvec] = gmres(Ax,b,50,1e-6,3,M_bs_y);

  %  compute new guess
  X = X + delta_X;
  abs_err = norm(delta_X)


  %  compute residual

  unew = X(1:end-2);
  omega = X(end-1)
  t_end = X(end)

  dt

  tspan = [0 t_end];
  Nstep=round(t_end/dt);
  dt_now=t_end/Nstep;


  [phi, p_full] =  TS_3x_trunc(tspan,unew,Ra,dt_now,K,Nc);
  gamma = rotation_trunc(unew,omega,tspan,Nc,K);

  b = -[(phi-gamma) ; ...
          u_theta_0.'*(unew); ...
             u0.'*unew - u0.'*u0];

  New_res = norm(b,'inf')


    fprintf('\n \n')
    
  New_iter = New_iter + 1;
  
  test =((New_res>abs_tol) & (New_iter < imax));  
  
  disp(['parameter = ' num2str(Ra) '        ' ' Newton iteration = ' num2str(New_iter)])
  disp([' residual = ' num2str(New_res),' abs error = ' num2str(abs_err) ])
  disp(['flag = ' int2str(flag_gmres),'    ' ' gmres_iter = ' int2str(gmres_iter),'    ' ' residual = ' num2str(relres)])

    New_res_vec(New_iter) = New_res;
    New_abs_vec(New_iter) = abs_err;
    Gmres_rel_vec(New_iter) = relres;
    Gmres_flag_vec(New_iter) = flag_gmres;

end
tic

fprintf('\n  ')

u_star = X(1:end-2);
nrm_u_st = norm(u_star);

phi_period = phi;
period = t_end;

format long


 %%%%%%%  eigenvalue computation %%%%%%%%

gmit=0;
fprintf('\n eig iteration:  ')

opts.disp=2

[V,D,flag_eigs] = eigs(@(v) matvec_prod_vac_eig(u_star,phi_period,omega,t_end,dt,Ra,K,Nc,v),2*(Nc+1)*(2*K+1),n_eigs,'lm',opts);    

mu = diag(D);
fprintf('\n')
disp(['eigs flag = ' , int2str(flag_eigs)])

fprintf('\n eigs: \n')

disp([mu abs(mu)])
eig_time = toc;

      
