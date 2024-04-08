function [gamma] = rotation_trunc(u,omega,tspan,Nc,K)

kwave = 0:K;
tau = diff(tspan);
A = zeros(Nc+1,length(kwave));
for k = 1:Nc+1
    A(k,:) = 1i*kwave*omega*tau;
end
B = exp(-A);

[psi2,phi] = vec_2_mat_trunc(u,Nc,K);
psi2m = my_fft(psi2,2*K+1);
phim  = my_fft(phi,2*K+1);

psi2m_rot = B.*psi2m;
phim_rot  = B.*phim;

psi2_now = my_ifft(psi2m_rot,2*K+1);
phi_now  = my_ifft(phim_rot,2*K+1);

gamma = mat_2_vec_trunc(psi2_now,phi_now,Nc,K);
end



