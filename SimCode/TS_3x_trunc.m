function [U_t, U_perturbation] = TS_3x_trunc(tspan,Up_t,Ra,dt,K,Nc)

load Tq_vs_psi2_Nc24_K32_rr_0p56_k_decay_value_50.mat
%load('base_state_charge_rr_p56Nc24.mat')
load('base_state_charge_Nc24_rr_0p56.mat')

%load Tq_vs_psi2_Nc48_K64_rr_0p56_k_decay_value_50.mat
%load base_state_charge_Nc48_rr_0p56.mat

%dt  = 1.0e-04;   % size of the time step:
Pr  = 75.8;       % Prandtl number:
Re  = 0.231;      % Reynolds number:
rr  = 0.56;       % relative ratio
%K   = 32;        % number of Fourier modes in the angular direction
%Nc  = 24;        % number of meshpoints in the radial direction:


%if nargin < 4,    dt  = 1.0e-04;    end  % size of the time step:
%if nargin < 5,    Pr  = 75.8;       end  % Prandtl number: 
%if nargin < 6,    Re  = 0.231;      end  % Reynolds number:
%if nargin < 7,    rr  = 0.56;        end  % relative ratio
%if nargin < 8,    K   = 32;         end  % number of Fourier modes in the angular direction
%if nargin < 9,    Nc  = 24;         end  % number of meshpoints in the radial direction:

Omega = Re*Pr*(1-rr)/rr;    % the rotation rate of the inner electrode.  Need this to compute the base state.
Nrun = round((tspan(2)-tspan(1))/dt);
% The following are used by the script diagnostics.m which computes diagnostic quantities such 
% as the kinetic energy density, the enstrophy density, and the Nusselt number all as a function of time.
% charge relaxation timescale
eps0        = 8.85*10^(-12);    %permittivity of free space
sigma3      = 1.3*10^(-7);      % from page 3624 of Daya Phys Fluids 11(1999)
one_layer_s = 3.6*10^(-9);      % each layer is 3.6nm
num_layers  = 75;               % from Stephen
s           = num_layers*one_layer_s;
sigma       = sigma3*s;
r_out       = 6.4*10^(-3);      % from page 3624 of Daya Phys Fluids 11(1999)
r_in        = rr*r_out;         % by definition of the radius ratio rr
d           = r_out - r_in;     % by definition of d
tau_q       = eps0*d/sigma;     % the charge relaxation time

% In the time-stepping, we need to update the vorticity.  This uses BCs4dr_phi
BCs4dr_phi  = [0;0];
% to turn off the nonlinear terms for diagnostic reasons, set Nonlinear to zero.
Nonlinear   = 1;

Nprint=1;
phi_BC=1;


xfact=5;
omags=6;
%xfact=2;
%omags=2;



% ==============================================================
% Create Grids and the differentiation Matrices
% ==============================================================
[T_inv,c_int,D,D2,yi,I,It,theta,riprime] = make_grids(Nc,K,rr);
% ==============================================================
%   initialize various variables.  This creates a bunch of zero arrays
% ==============================================================
[Ur,Utheta,Phi,psi2,q,W,Flux,time,E_kin,E_enstr,Nu,psi2_norm,q_norm,phi_norm,w_norm] = initialize(Nc,K,Nprint);
% =============================================================
%  Here's the T-matrix
% =============================================================
%k_max=50;
%filename = strcat('Tq_vs_psi2_Nc',num2str(Nc),'_K',num2str(K),'_rr_',num2str(rr),'_k_decay_value_',num2str(k_max),'.mat');
% If you want to compute it using the symbolic toolbox then use the function
% qm_vs_Psi2m.m .  But be warned --- this is slow!!
% the transformation function between hat{q_m} and hat{psi2_m}
%Nk=50000;
% Tq_vs_psi2 = qm_vs_Psi2m(K, Nc,Nk,k_max, rr);
% the spectral integration in the radial direction...
% save(filename,'Tq_vs_psi2');
% or load the precomputed (Nc+1)x(Nc+1)x(K+1) matrix Tq_vs_psi2
%load(filename)
%
% =========================================================
% Compute the base state and its derivatives (needed for the Jacobians)
% =========================================================
%filename = strcat('base_state_charge_Nc',num2str(Nc),'_rr_',num2str(rr),'.mat');
% compute the base state charge and potential using the symbolic toolbox --- this is a little slow!
% [Psi20,Psi20_r,Q0,Q0_r] = base_state_charge(rr, riprime);
% save(filename,'Psi20','Psi20_r','Q0','Q0_r')
% or load them directly from a file
%load(filename)
% load base_state_charge_Nc24_rr_0.56.mat

% As given, the base state quantities are independent of theta; they're
% (Nc+1) x 1 matrices representing the radial component.  We want to create
% (Nc+1) x (2K+1) matrices which will represent the base state quantities
% in the annulus for use in reconstructing the unknowns in the annulus
% (find_fields.m) as well as in computing diagnostic quantities (diagnostics.m
% and norms.m)
Psi20 = Psi20*ones(1,2*K+1);
Psi20_r = Psi20_r*ones(1,2*K+1);
Q0 = Q0*ones(1,2*K+1);

% compute the base state stream function and vorticity so they can be
% compared to the computed solution as needed.  If the inner electrode is
% moving in a time-dependent manner then we'll want to compute this at
% every timestep.
[Phi0,Phi0_r,W0] = base_state_flow(Omega, rr, riprime);
Phi0 = Phi0*ones(1,2*K+1);
Utheta0 = -Phi0_r*ones(1,2*K+1);
W0 = W0*ones(1,2*K+1);
Flux0 = -Phi0_r;


% =======================================================
%   Initial Conditions and start-up
% =======================================================
% initialize arrays that would otherwise dynamically grow:
% [Ur_phy,Utheta_phy,phi_phy,psi2_phy,q_phy,w_phy,time,E_kin,E_enstr,Nu,psi2_norm,q_norm,phi_norm,w_norm] = initialize(Nc,K,Nrun);
% for random initial data:
% random_initial_data_PT
%[psi2_now,psi2m_now,q_now,qm_now,w_now,wm_now,phi_now,phim_now] = random_initial_data(Tq_vs_psi2,D,D2,riprime,theta,Nc,K,InitBCs_amp);

%[psi2m_now,qm_now,wm_now,phim_now] =  reshape_IC_time_stepper(Up_in,Nc,K);
%[psi2_now,q_now,w_now,phi_now] =  vec_2_mat(Up_in,Nc,K);

% set up initial conditions

[psi2m_now,qm_now,wm_now,phim_now,psi2_now,q_now,w_now,phi_now] = truncvec_2_allmat(Up_t,Tq_vs_psi2,D2,D,Omega,rr,riprime,Nc,K);


 Psi2_now=psi2_now+Psi20;
 Q_now=q_now+Q0;
 W_now=w_now+W0;
 Phi_now=phi_now+Phi0;

 Wm_now = my_fft(W_now,2*K+1);      
 Phim_now = my_fft(Phi_now,2*K+1); 


%[psi2_now,q_now,w_now,phi_now] =  vec_2_mat(U0,Nc,K);
% if you want to use the same initial data for multiple runs then load it
% from a file:
% load rand_ID.mat
% t_start = 0;
% If you want to continue a run from data that you saved already in a file
% [Ra,Pr,Re,rr,psi2_now,psi2m_now,q_now,qm_now,phi_now,phim_now,w_now,wm_now,t_start] = load_ID('Pr10_Ra70_Re0p2_rr0p56.mat');

% Why am I using my own fft and ifft you may well ask?  Well the arrays in
% this code are being stored in the opposite way that matlab's fft/ifft are
% expecting.  And so I have to do transposes and stuff like that.  This is
% a vestige of how Peichun wrote the code; I think she stored the data in
% this particular way so that she could use polar3d.m to plot the data.

% in case you want to see the spectra of the initial data
% plot_spectra
% pause(1)


dt = dt/(xfact^omags);


% set the inner rotation speed to Omega*r_in
in_speed = Omega*rr/(1-rr);
delta=dt*Pr;

% precompute the influence matrix, needed for the vorticity/fluid potential
% solve.
[wm1,phim1,wm3,phim3,Mi] = influence_matrix(D,D2,riprime,K,delta,It);

% fprintf('//Ra = %.2f:\n', Ra);
% fprintf('//Pr = %.2f:\n', Pr);
% fprintf('//Re = %.2f:\n', Re);

% we have Utheta(r,theta,t).  We want to compute its average in theta:
%   flux(r,t) = 1/(2pi) \int_0^2pi Utheta(r,theta,t) dtheta.  This is done
% with spectral accuracy using the trapezoid rule.  dtheta = 2*pi/(2*K+1)
%   flux(r,t) = 1/(2pi) sum(Utheta) (2pi)/(2K+1) = sum(Utheta)/(2K+1).
% Recall that Utheta starts at theta=0 and goes up to theta=2*pi-dtheta. The
% trapezoid rule which would usually have Utheta(0)/2 + ... + Utheta(2pi)/2
% is replaced by Utheta(0)+ ...+ Utheta(2pi-dr) .  Also, I have to be careful
% to get the summation in the right direction (columns versus rows) and all
% that stuff.  The sum command lets me say whether I want to sum over the
% first dimension (rows) or over the second dimension (columns).

drPhim = 2*D*Phim_now;
drPhi = my_ifft(drPhim,2*K+1);
Utheta = - drPhi;

flux_now = sum(Utheta(:,:),2)/(2*K+1);
%Flux(:)=flux_now;

% ===================== Begin the time iteration ============================
% we're using a three-level timestepping scheme.  So we start with one
% step of a first-order accurate two-level scheme.  And then continue on
% with the three-level scheme.


% need to initialize these variables.
psi2m_new = zeros(Nc+1,K+1);
qm_new = zeros(Nc+1,K+1);
Phim_new = zeros(Nc+1,K+1);
Wm_new = zeros(Nc+1,K+1);


for jjfirst=1:1  % initial steps with first order method

% compute the initial terms that go into the nonlinearity of the PDEs for
% psi2, q, phi, and w
[Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(psi2m_now,...
    Wm_now,Phim_now,qm_now,Q0_r,Psi20_r,riprime,Nonlinear,D,Nc,2*K);
% first we define the nonlinear terms that're on the RHS of the flux PDE
flux_nonl_now = flux_nonlinearity(Pr,Ra,psi2m_now,qm_now,Q0,Phim_now,riprime,D,K);


if (jjfirst==1)
  psi2m_old = psi2m_now;
  qm_old = qm_now;
  Wm_old = Wm_now;
  Phim_old = Phim_now;
  flux_old = flux_now;
  flux_nonl_old = flux_nonl_now;
  Jqphi_hat_old = Jqphi_hat_now;
  Jwphi_hat_old = Jwphi_hat_now;
  Jpsi2q_hat_old = Jpsi2q_hat_now;

  psi2m_first = psi2m_now;
  qm_first = qm_now;
  Wm_first = Wm_now;
  Phim_first = Phim_now;
  flux_first = flux_now;
  flux_nonl_first = flux_nonl_now;
  Jqphi_hat_first = Jqphi_hat_now;
  Jwphi_hat_first = Jwphi_hat_now;
  Jpsi2q_hat_first = Jpsi2q_hat_now;

end


% take one step using first-order time-stepping
% precompute the influence matrix, needed for the vorticity/fluid potential
% solve.
delta = dt*Pr;
[wm1,phim1,wm3,phim3,Mi] = influence_matrix(D,D2,riprime,K,delta,It);
[psi2m_new,qm_new,Phim_new,Wm_new,flux_new] = first_order_time_step(psi2m_new,...
    qm_now,qm_new,Phim_new,Wm_now,Wm_new,flux_now,flux_nonl_now,Jqphi_hat_now,...
    Jwphi_hat_now,Jpsi2q_hat_now,D,D2,K,riprime,...
    dt,Tq_vs_psi2,Nc,T_inv,c_int,wm1,wm3,phim1,phim3,Mi,phi_BC,Ra,Pr,in_speed,It);
% update so we can take the next timestep.

psi2m_now = psi2m_new;
qm_now = qm_new;
Wm_now = Wm_new;
Phim_now = Phim_new;

end


for jjlps=1:omags

% compute the Jacobian terms pseudospectrally
[Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(psi2m_now,...
    Wm_now,Phim_now,qm_now,Q0_r,Psi20_r,riprime,Nonlinear,D,Nc,2*K);
flux_now = flux_new;

% compute the flux nonlinearity
flux_nonl_now = flux_nonlinearity(Pr,Ra,psi2m_now,qm_now,Q0,Phim_now,...
    riprime,D,K);

% precompute the influence matrix, needed for the vorticity/fluid potential
% solve.  We have to do this again because for the first-order scheme we
% were trying to invert a different matrix than for the second-order scheme
delta = (2/3)*dt*Pr;
[wm1,phim1,wm3,phim3,Mi] = influence_matrix(D,D2,riprime,K,delta,It);


for jjsec=2:xfact  % initial steps with 2nd order method

    [psi2m_new,qm_new,Phim_new,Wm_new,flux_new] = second_order_time_step(psi2m_new,...
        qm_old,qm_now,qm_new,Phim_new,Wm_old,Wm_now,Wm_new,flux_now,flux_old,flux_nonl_now,flux_nonl_old,...
        Jqphi_hat_old,Jqphi_hat_now,Jwphi_hat_old,Jwphi_hat_now,Jpsi2q_hat_old,...
        Jpsi2q_hat_now,D,D2,K,riprime,dt,Tq_vs_psi2,Nc,T_inv,c_int,wm1,wm3,phim1,phim3,...
        Mi,phi_BC,Ra,Pr,in_speed,It);
    % update so we can take the next timestep.
    psi2m_old = psi2m_now;
    qm_old = qm_now;
    Wm_old = Wm_now;
    Phim_old = Phim_now;
    flux_old = flux_now;
    flux_nonl_old = flux_nonl_now;
    Jqphi_hat_old = Jqphi_hat_now;
    Jwphi_hat_old = Jwphi_hat_now;
    Jpsi2q_hat_old = Jpsi2q_hat_now;
    psi2m_now = psi2m_new;
    qm_now = qm_new;
    Wm_now = Wm_new;
    Phim_now = Phim_new;
    % compute the Jacobian terms pseudospectrally
    [Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(...
        psi2m_now,Wm_now,Phim_now,qm_now,Q0_r,Psi20_r,riprime,Nonlinear,D,...
        Nc,2*K);
    flux_now = flux_new;
    % compute the flux nonlinearity
    flux_nonl_now = flux_nonlinearity(Pr,Ra,psi2m_now,qm_now,Q0,Phim_now,...
        riprime,D,K);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%  now have 1 step at dt (computed from xfact steps at dt/xfact)
% %%  now go to next higher order of magnitude

dt = dt*xfact;

  psi2m_old = psi2m_first;
  qm_old = qm_first;
  Wm_old = Wm_first;
  Phim_old = Phim_first;
  flux_old = flux_first;
  flux_nonl_old = flux_nonl_first;
  Jqphi_hat_old = Jqphi_hat_first;
  Jwphi_hat_old = Jwphi_hat_first;
  Jpsi2q_hat_old = Jpsi2q_hat_first;

end

% compute the Jacobian terms pseudospectrally
[Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(psi2m_now,...
    Wm_now,Phim_now,qm_now,Q0_r,Psi20_r,riprime,Nonlinear,D,Nc,2*K);
flux_now = flux_new;
% compute the flux nonlinearity
flux_nonl_now = flux_nonlinearity(Pr,Ra,psi2m_now,qm_now,Q0,Phim_now,...
    riprime,D,K);

% precompute the influence matrix, needed for the vorticity/fluid potential
% solve.  We have to do this again because for the first-order scheme we
% were trying to invert a different matrix than for the second-order scheme
delta = (2/3)*dt*Pr;
[wm1,phim1,wm3,phim3,Mi] = influence_matrix(D,D2,riprime,K,delta,It);



for jj=2:Nrun
    [psi2m_new,qm_new,Phim_new,Wm_new,flux_new] = second_order_time_step(psi2m_new,...
        qm_old,qm_now,qm_new,Phim_new,Wm_old,Wm_now,Wm_new,flux_now,flux_old,flux_nonl_now,flux_nonl_old,...
        Jqphi_hat_old,Jqphi_hat_now,Jwphi_hat_old,Jwphi_hat_now,Jpsi2q_hat_old,...
        Jpsi2q_hat_now,D,D2,K,riprime,dt,Tq_vs_psi2,Nc,T_inv,c_int,wm1,wm3,phim1,phim3,...
        Mi,phi_BC,Ra,Pr,in_speed,It);
    % update so we can take the next timestep.
    psi2m_old = psi2m_now;
    qm_old = qm_now;
    Wm_old = Wm_now;
    Phim_old = Phim_now;
    flux_old = flux_now;
    flux_nonl_old = flux_nonl_now;
    Jqphi_hat_old = Jqphi_hat_now;
    Jwphi_hat_old = Jwphi_hat_now;
    Jpsi2q_hat_old = Jpsi2q_hat_now;
    psi2m_now = psi2m_new;
    qm_now = qm_new;
    Wm_now = Wm_new;
    Phim_now = Phim_new;
    % compute the Jacobian terms pseudospectrally
    [Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(...
        psi2m_now,Wm_now,Phim_now,qm_now,Q0_r,Psi20_r,riprime,Nonlinear,D,...
        Nc,2*K);
    flux_now = flux_new;
    % compute the flux nonlinearity
    flux_nonl_now = flux_nonlinearity(Pr,Ra,psi2m_now,qm_now,Q0,Phim_now,...
        riprime,D,K);

end


[Psi2,psi2,Psi2_r,Q,q,Phi,phi,W,w,Ur,Utheta] ...
    = find_fields(Psi20,Psi20_r,Q0,psi2m_now,qm_now,Phi0,Phim_now,W0,Wm_now,riprime,D,K);


U_t = mat_2_vec_trunc(psi2,phi,Nc,K);
U_perturbation = mat_2_vec(psi2,q,w,phi,Nc,K);


[psi2m_now,qm_now,wm_now,phim_now,psi2_now,q_now,w_now,phi_now] = truncvec_2_allmat(U_t,Tq_vs_psi2,D2,D,Omega,rr,riprime,Nc,K);


end


    % in case you want to see what the vorticity looks like right now
    %     figure(2)
    %     polar3d(W(:,:,ii+1),0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
    %     figure(2)

    %     % in case you want to see the power spectrum
    %     plot_spectra(T_inv,psi2_now,K)
    %     pause(1)
    %     plot_spectra(T_inv,q_now,K)
    %     pause(1)
    %     plot_spectra(T_inv,Phi_now,K)
    %     pause(1)
    %     plot_spectra(T_inv,W_now,K)
    %     pause(1)

    %     % in case you want to see all four uknowns right now
    %     figure(1)
    %     plot_soln
    %     figure(1)
    %
    %     figure(2)
    %     polar3d(psi2(:,:,ii+1),0,2*pi,rr,1,1,'contour');
    %     title('psi2')
    %     figure(2)



