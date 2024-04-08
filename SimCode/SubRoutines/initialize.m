function [Ur,Utheta,Phi,psi2,q,W,Flux,time,E_kin,E_enstr,Nu,psi2_norm,q_norm,phi_norm,w_norm] = initialize(Nc,K,Nprint);

E_enstr = zeros(1,Nprint+1);
E_kin = zeros(1,Nprint+1);
phi_norm = zeros(1,Nprint+1);
psi2_norm = zeros(1,Nprint+1);
q_norm = zeros(1,Nprint+1);
time = zeros(1,Nprint+1);
w_norm = zeros(1,Nprint+1);

Flux = zeros(Nc+1,Nprint+1);
Nu = zeros(Nc+1,Nprint+1);

Phi = zeros(Nc+1,2*K+1,Nprint+1);
Ur = zeros(Nc+1,2*K+1,Nprint+1);
Utheta = zeros(Nc+1,2*K+1,Nprint+1);
W = zeros(Nc+1,2*K+1,Nprint+1);
psi2 = zeros(Nc+1,2*K+1,Nprint+1);
q = zeros(Nc+1,2*K+1,Nprint+1);
             
  
