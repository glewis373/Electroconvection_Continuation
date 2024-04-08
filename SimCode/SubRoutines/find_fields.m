function [Psi2,psi2,Psi2_r,Q,q,Phi,phi,W,w,Ur,Utheta] ...
    = find_fields(Psi20,Psi20_r,Q0,psi2m,qm,Phi0,Phim,W0,Wm,riprime,D,K)
% Given the base states (indicated with "0"s in their names) 
% the code uses the deviations in spectral space (indicated by lower case) 
% and produces the full solution in physical space by adding them together
% (done for psi2 and q)
%
% Also, the code takes the full solution in spectral space and produces 
% the deviations by subtracting off the base state (done for Phi and W).
%
% In addition, the radial and angular velocities are computed.

% r-component of the velocity is 1/r dPhi/dtheta.
testm = Phim*diag(1i*(0:1:K));
dPhi = my_ifft(testm,2*K+1);
Ur = diag(1./riprime)*dPhi;
% theta-component of the velocity is -dPhi/dr.  Factor of 2 is from the
% chain rule.
drPhim = 2*D*Phim;
drPhi = my_ifft(drPhim,2*K+1);
Utheta = - drPhi;


% construct the other quanties of interest
% the deviation of the potential from base state in physical space
psi2 = my_ifft(psi2m,2*K+1);  
% the deviation of the charge from base state in physical space
q = my_ifft(qm,2*K+1); 
% Now add these deviations to the base state to get the potential and charge
% in physical space.
Psi2 = psi2 + Psi20;
Q = q + Q0;

% the deviation of the stream function from base state in physical space
Phi = my_ifft(Phim,2*K+1); 
% the deviation of the vorticity from base state in physical space
W = my_ifft(Wm,2*K+1);
% Now subtract off the base state to get the deviations from the base state 
% of the stream function and vorticity in physical space.
phi = Phi-Phi0;
w = W-W0;

% some things needed in order to compute the Nusselt number:
drpsi2m = 2*D*psi2m;
drpsi2 = my_ifft(drpsi2m,2*K+1);
Psi2_r = drpsi2 + Psi20_r;
