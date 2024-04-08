function [psi2,psi2m,q,qm,W,Wm,Phi,Phim] = base_state_initial_data(Omega,rr,riprime,Nc,K)

% Take the potential and charge distribution as zero
psi2 = zeros(Nc+1,2*K+1);
psi2m = zeros(Nc+1,K+1);
q = zeros(Nc+1,2*K+1);
qm = zeros(Nc+1,K+1);

% Take the stream function and vorticity base state
% compute the base state stream function and vorticity
[Phi,Phi_r,W] = base_state_flow(Omega, rr, riprime); 
Phi = Phi*ones(1,2*K+1);
W = W*ones(1,2*K+1);
Phim = my_fft(Phi,2*K+1);
Wm = my_fft(W,2*K+1);
