function [Ax] = matvec_prod_rel_equi(unew,phi,u0,u_theta_0,omega,tspan,dt,Ra,K,Nc,Dv)

global gmit

%K   = 32;         end  % number of Fourier modes in the angular direction
%Nc  = 24;         end  % number of meshpoints in the radial direction:

alpha = 0;
epsilon = 1e-4;
Dv1 = Dv(1:end-2);
Dv2 = Dv(end-1);
Dv3 = Dv(end);


gmit=gmit+1;
fprintf('%d;  ',gmit)

%nrm_Dv1=norm(Dv1)


% forward finite difference for the Jacobian Matrix


tau=tspan(2)+epsilon*Dv3;
tspan_e=[0 tau];

Nstep=round(tau/dt);
dt_now=tau/Nstep;


vin  = unew +epsilon*Dv1;
[phi_e_total, p_full] = TS_3x_trunc(tspan_e,vin,Ra,dt_now,K,Nc);
phi_e = phi_e_total(:,end);
J = (phi_e-phi)/epsilon;


gamma_0 = rotation_trunc(unew,omega,tspan,Nc,K);
gamma_pert = rotation_trunc(vin,omega+epsilon*Dv2,tspan_e,Nc,K);

d_rot = (gamma_pert - gamma_0)/epsilon;


Ax1 = J - d_rot;
Ax2 = u_theta_0.'*Dv1 + alpha*Dv2;
Ax3 = u0.'*Dv1;
Ax =[Ax1;Ax2;Ax3];
end
