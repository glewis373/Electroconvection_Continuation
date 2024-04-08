function [Wmnew,Phimnew] = first_order_vorticity_update(Wmnew,Phimnew,Wm_now,flux_new,D,wm1,wm3,phim1,phim3,Mi,phi_BC,Ra,Pr,in_speed,dt,It,Lop,m,Nc,T_inv,c_int,K,Jwphi_hat_now,Jpsi2q_hat_now)

% One step scheme, doing the nonlinear term explicitly and the diffusion
% term implicitly
%
% (W_new-W_now)/dt + Jwphi_now = Pr L W_new + Pr Ra Jpsi2q_now
%
% W_new-W_now + dt Jwphi_now = dt Pr L W_new + dtPr Ra Jpsi2q_now
%
% W_new - dt Pr L W_new = W_now - dt Jwphi_now  + dt Pr Ra Jpsi2q_now
%
% (I - dt Pr L ) W_new = W_now - dt Jwphi_now  + dt Pr Ra Jpsi2q_now

% See second_order_vorticity_update.m for description of how method works
if m == 0 
    % have boundary conditions on the zero mode
    % here is the (negative) of the velocity at the inner radius
    h1 = - in_speed;
    if phi_BC == 1
        % here is the total flux
        % g1 = Chebyshev_Int(Nc,flux_new)/2;
        g1 = dot(c_int,T_inv*flux_new)/2;
    else
        g1 = 0;
    end
    % because the boundary conditions are imposed in Fourier space, not in
    % real space, we need to multiply by the factor the the FFT brings to
    % the party.
    h1 = (2*K+1)*h1;
    g1 = (2*K+1)*g1;
else
    h1 = 0;
    g1 = 0;
end
g2 = 0;
h2 = 0;

% solve the "step 1 problem"
A_full = (It - dt*Pr*Lop(:,:));
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on the vorticity
Wm2(1,m+1)=0;
Wm2(Nc+1,m+1)=0;
% We use these to adjust the RHS...
B = Wm_now(2:Nc, m+1) - dt*Jwphi_hat_now(2:Nc, m+1)...
    + dt*Pr*Ra*Jpsi2q_hat_now(2:Nc, m+1);
RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[Wm2(1,m+1);Wm2(Nc+1,m+1)];
Wm2(2:Nc,m+1) = A_trunc\RHS;
% Now compute the potential from wm2...
A_full = Lop;
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on the potential
Phim2(1,m+1)=g2;
Phim2(Nc+1,m+1)=g1; 
B = -Wm2(2:Nc,m+1);
RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[Phim2(1,m+1);Phim2(Nc+1,m+1)];
Phim2(2:Nc,m+1) = A_trunc\RHS;

% compute the coefficients needed to solve the "step 2 problem"
BCs4w (:, m+1) = Mi(:,:,m+1)\( [h1;h2] - [2*D(Nc+1, :)*Phim2(:,m+1); 2*D(1, :)*Phim2(:, m+1)] );
% define the solution as the sum of the step 1 and step 2 problems
%Wmnew(:,m+1)=Wm2(:,m+1)+BCs4w(1,m+1)*wm1(:,m+1)+BCs4w(2,m+1)*wm3(:,m+1);
Phimnew(:,m+1)=Phim2(:,m+1)+BCs4w(1,m+1)*phim1(:,m+1)+BCs4w(2,m+1)*phim3(:,m+1);

Wmnew(:,m+1)=-Lop*Phimnew(:,m+1);

