function flux_new = second_order_flux_update(flux_now,flux_old,flux_nonl_now,flux_nonl_old,riprime,in_speed,Pr,It,D,D2,dt,Nc)

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
% the deviation.  (A moot point since we've set Phi_0 to zero.)
% That is, u_r := 1/r dPhi/dtheta and u_theta := - dPhi/dr.

% The two-step time-stepping is explicit on the nonlinear terms and implicit on the
% linear terms:
%
% (3*flux_new - 4*flux_now + flux_old)/(2*dt) 
%             = Pr Delta_r flux_new - Pr/r^2 flux_new + 2*flux_nonl_now - flux_nonl_old 
%
% where flux_nonl := - Pr Ra/r mean(Q dPsi2/dtheta)
%                      - mean( u_r d/dtheta u_theta) - 1/r mean( u_r u_theta)
% all at the current time-step; it's computed using flux_nonlinearity.m
%
% So we have 
%   3*flux_new - 4*flux_now + flux_old = 2*dt*Pr*Delta_r flux_new 
%                -2*dt*Pr/r^2 flux_new + 2*dt*(2*flux_nonl_now - flux_nonl_old)
% and so
%   flux_new - 4/3*flux_now + 1/3 flux_old = 2/3*dt*Pr*Delta_r flux_new
%                   -2/3*dt*Pr/r^2 flux_new + 2/3*dt*(2*flux_nonl_now - flux_nonl_old)
% and so
%   flux_new - 2/3*dt*Pr*Delta_r flux_new + 2/3*dt*Pr/r^2 flux_new 
%              = 4/3*flux_now - 1/3 flux_old + 2/3*dt*(2*flux_nonl_now - flux_nonl_old)
%
% A flux_new = B
%
% where A := I - 2*dt/3*Pr*(Lop - 1/r^2*I)    
% and B = 4/3*flux_now -1/3*flux_old + 2/3*dt*(2*flux_nonl_now-flux_nonl_old)

A_full = 4*D2 + 2*diag(1./riprime)*D - diag(1./(riprime.^2));
A_full = (2/3)*dt*Pr*A_full;
A_full = It - A_full;
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on flux_new
flux_new(1,1)=0; % outer electrode is at rest
flux_new(Nc+1,1)=in_speed; % inner electrode has speed in_speed
% We use these to adjust the RHS...
B = (4/3)*flux_now(2:Nc)-(1/3)*flux_old(2:Nc) + (2/3)*dt*(2*flux_nonl_now(2:Nc)-flux_nonl_old(2:Nc));
RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[flux_new(1,1);flux_new(Nc+1,1)];
flux_new(2:Nc,1) = A_trunc\RHS;
