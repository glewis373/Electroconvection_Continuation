function [Av] = matvec_prod_vac_eig(u_star,phi_period,omega,period,dt,Ra,K,Nc,v)

global gmit

%K   = 32;         end  % number of Fourier modes in the angular direction
%Nc  = 24;         end  % number of meshpoints in the radial direction:

epsilon = 1e-4;

gmit=gmit+1;
fprintf('%d;  ',gmit)

Nstep=round(period/dt);
dt_now=period/Nstep;
tspan=[0 period];

vin    = u_star + epsilon*v;
phi_ue = TS_3x_trunc([0, period],vin,Ra,dt_now,K,Nc);

gamma_0 = rotation_trunc(phi_period,-omega,tspan,Nc,K);
gamma_pert = rotation_trunc(phi_ue,-omega,tspan,Nc,K);

Av=(gamma_pert - gamma_0)/epsilon;



end
