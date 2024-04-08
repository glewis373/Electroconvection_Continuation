function flux_nonl= flux_nonlinearity(Pr,Ra,psi2m,qm,Q0_phy,Phim,riprime,D,K) 
% We have the PDE
%
%  d/dt flux = Pr Delta_r flux - Pr/r^2 flux - Pr Ra/r mean(Q dPsi2/dtheta)
%                  - mean( u_r d/dtheta u_theta) - 1/r mean( u_r u_theta)
%
% where flux(r,t) := 1/2pi \int u_theta(r,theta,t) dtheta
% and   Delta_r := the radial Laplacian = d^2/dr^2 + 1/r d/dr
% and   mean(f(r,theta,t)) := 1/2pi \int f(r,theta,t) dtheta
% and   Q and Psi2 are the physical charge distribution, not the deviation
% and   u_r and u_theta are the physical components of the velocity, not
% the deviation.  
% That is, u_r := 1/r dPhi/dtheta and u_theta := - dPhi/dr.
%
% In this subroutine, I want to compute the terms that will show up on the
% RHS.  For both the first-order and second-order time-stepping I'll handle
% the Pr Delta_r flux and - Pr/r^2 flux terms implicitly.  And I'll need
%
%   flux_nonl := - Pr Ra/r mean(Q dPsi2/dtheta)
%                  - mean( u_r d/dtheta u_theta) - 1/r mean( u_r u_theta)

% r-component of the velocity is 1/r dPhi/dtheta.
testm = Phim*diag(1i*(0:1:K));
dPhi = my_ifft(testm,2*K+1);
u_r = diag(1./riprime)*dPhi;
% theta-component of the velocity is -dPhi/dr.  Factor of 2 is from the
% chain rule.
drPhim = 2*D*Phim;
drPhi = my_ifft(drPhim,2*K+1);
u_theta = - drPhi;

% Start constructing the nonlinear term 
% step 1: build
%   flux_nonl = - Pr Ra/r Q dPsi2/dtheta - u_r d/dtheta u_theta - 1/r u_r u_theta
% step 2: average out the angular stuff:
%   flux_nonl = mean(flux_nonl)

% first the -1/r u_r u_theta term
flux_nonl = - diag(1./riprime)*(u_r.*u_theta);

% now the - u_r d/dtheta u_theta term
% Need d/dtheta u_theta = d/dtheta (- d/dr Phi) 
testm = -drPhim*diag(1i*(0:1:K));
dthetadru_theta = my_ifft(testm,2*K+1);
flux_nonl = flux_nonl - u_r.*dthetadru_theta;

% now the -Pr Ra/r Q dPsi2/dtheta term.  Note that we really have 
%      dPsi2/dtheta = d/dtheta (psi2_phy + Psi20_phy)
%                   = d/dtheta (psi2_phy)
% and so we only need to compute the theta-derivative of the deviation
% psi2.  
testm = psi2m*diag(1i*(0:1:K));
dthetaPsi2 = my_ifft(testm,2*K+1);
% the deviation of the charge from base state in physical space
q_phy = my_ifft(qm,2*K+1); 
% Now add this deviation to the base state to get the charge in physical
% space.
Q_phy = q_phy + Q0_phy;
flux_nonl = flux_nonl - Pr*Ra*diag(1./riprime)*(Q_phy.*dthetaPsi2);

% now average in theta:
flux_nonl = sum(flux_nonl,2)/(2*K+1);


