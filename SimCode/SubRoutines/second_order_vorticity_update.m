function [Wm_new,Phim_new] = second_order_vorticity_update(Wm_new,Phim_new,Wm_now,Wm_old,flux_new,D,wm1,wm3,phim1,phim3,Mi,phi_BC,Ra,Pr,in_speed,dt,It,Lop,m,Nc,T_inv,c_int,K,Jwphi_hat_now,Jwphi_hat_old,Jpsi2q_hat_now,Jpsi2q_hat_old)
% We're trying to solve
% 
%           (I-delta L) w = RHS
%                  L  Phi = -w
% 
% subject to the boundary conditions 
%       Phi(-1)=g1, Phi(1)=g2, Phi'(-1)=h1, Phi'(1)=h2
% w/ g2=h2=0 and h1=-rin*Omega
% if m = 0 then
%      g1 = \int_rin^rout u_theta(r,theta0,t) dr = Phi(-1)-Phi(1)
% if m \neq 0 then
%      g1 = 0
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
%
%
% Because this has 4 BCs on Phi and no BCs on w, we use the influence
% matrix method to handle the problem; this takes major advantage of the
% linearity of the problem.
%
% Note: the problem is of the form A1 w = RHS, A2 Phi = -w.  So first solve
% the A1 w = RHS problem to find w.  Now that w is known, we can solve the
% A2 Phi = -w problem for Phi.
%
% Step 1: solve the BVP
%          A1 w = RHS, A2 Phi = -w
% w(-1)=0, w(1)=0, Phi(-1)=g1, Phi(1)=g2.
% Call the solution of this BVP "w2" and "Phi2".
%
% Step 2: using the solution (w2,Phi2) define a new BVP
%          A1 w = 0, A2 Phi = -w
% Phi(-1)=0, Phi(1)=0, Phi'(-1)=h1-Phi2'(-1), Phi'(1)=h2-Phi2'(1).
%
% Phi := Phi-step1 + Phi-step2
% w := w-step1 + w-step2
%
% will solve the initial BVP.
%
% Note: the step 2 BVP is again a problem to solve --- all the BCs are on
% Phi.  To solve this, solve two auxiliary BVPS:
%
% step 2a: A1 w = 0, A2 Phi = -w with Phi(+/- 1) = 0, w(-1)=1, w(1) = 0
%        call the soln w1, phi1
% step 2b: A1 w = 0, A2 Phi = -w with Phi(+/- 1) = 0, w(-1)=0, w(1) = 1
%        call the soln w3, phi3
%
% By choosing the right linear combination of w1 and w3, of phi1 and phi3,
% we solve the BVP of step 2.  This involves solving M c = B where c is the
% coefficient vector.  NOTE: the matrix M depends on the solutions of steps
% 2a and steps 2b.  Which depend on dt and material properties but nothing
% else.  And so if dt is fixed then you can precompute the matrix M and
% just use c = M\B at each time step.  BUT if dt is not fixed in time, if
% it's an adaptive time-stepper, then either M needs to be computed at each
% time step or it needs to be precomputed for a fixed collection of
% timesteps.



% first step: solve
%           (I-delta L) w = RHS
%                    L  Phi = -w
% subject to the boundary conditions w(\pm 1)=0 and Phi(-1) = g1, Phi(1)=g2.
%
A_full = (It - 2/3*dt*Pr*Lop(:,:));
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on the vorticity
Wm2(1,m+1)=0;
Wm2(Nc+1,m+1)=0;
% We use these to adjust the RHS...
B = 4/3*Wm_now(2:Nc, m+1) - 1/3*Wm_old(2:Nc, m+1)...
    - 2/3*dt*(2*Jwphi_hat_now(2:Nc, m+1) - Jwphi_hat_old(2:Nc, m+1) )...
    + 2/3*dt*Pr*Ra*(2*Jpsi2q_hat_now(2:Nc, m+1) - Jpsi2q_hat_old(2:Nc, m+1)) ;
RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[Wm2(1,m+1);Wm2(Nc+1,m+1)];
Wm2(2:Nc,m+1) = A_trunc\RHS;
% Now compute the potential from Wm2...
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
%Wm_new(:,m+1)=Wm2(:,m+1)+BCs4w(1,m+1)*wm1(:,m+1)+BCs4w(2,m+1)*wm3(:,m+1);
Phim_new(:,m+1)=Phim2(:,m+1)+BCs4w(1,m+1)*phim1(:,m+1)+BCs4w(2,m+1)*phim3(:,m+1);


Wm_new(:,m+1)=-Lop*Phim_new(:,m+1);



