function [Uout] =  mat_2_vec_trunc(psi2_phy,phi_phy,Nc,K)
 
U1 = reshape(psi2_phy,(Nc+1)*(2*K+1),1);
U2 = reshape(phi_phy,(Nc+1)*(2*K+1),1);
Uout = [U1;U2];
end
