function [E_kin_dens,E_enstr_dens,Nu] ...
    = diagnostics(Ur,Utheta,W,Q,Psi20_r,Psi2_r,sigma,r,K,Nc)

% the area of the annulus
A = pi*(r(1)^2 - r(Nc+1)^2);
dtheta = 2*pi/(2*K+1);

% we want to approximate  \int_rin^rout \int_0^2pi f(r,theta) dtheta r dr

% integrate wrt angle using trapezoid rule to get 
%                F(r) = \int_0^2pi f(r,theta) dtheta
% the ",2" in the sum makes the summation go over the second index
f = 1/2*(Ur.^2 + Utheta.^2); % Kinetic energy
F = dtheta*sum(f,2);
% now we want to integrate r*F(r) wrt r using the Chebyshev integrator.
% The factor of 2 comes from change of coordinates, [r_in,r_out] -> [-1,1]
E_kin = Chebyshev_Int( Nc, r.*F )/2;
% divide by the area to get the kinetic energy density
E_kin_dens = E_kin/A;

f =  1/2*W.^2;  % Enstrophy = 1/2 vorticity^2.
F = dtheta*sum(f,2);
E_enstr = Chebyshev_Int( Nc, r.*F )/2;
E_enstr_dens = E_enstr/A;

% sigma3 = 1.3*10^(-7); % from page 3624 of Phys Fluids 11(1999)
% one_layer_s = 3.6*10^(-9); % each layer is 3.6nm
% num_layers = 75; % from Stephen
% s = num_layers*one_layer_s;
% sigma = sigma3*s;
arg1 = -sigma*Psi2_r + Q.*Ur;
arg2 = -sigma*Psi20_r;
Nu = sum(transpose(arg1))./sum(transpose(arg2)); 
% really the numerator and denominator should both be multiplied by dtheta 
% but then those would cancel out...


        
