function [psi2_now,phi_now] =  vec_2_mat_trunc(U_init,Nc,K)
len = (Nc+1)*(2*K+1);

psi2_now = reshape(U_init(1:len),(Nc+1),(2*K+1));
phi_now    = reshape(U_init(len+1:2*len),(Nc+1),(2*K+1));

end
