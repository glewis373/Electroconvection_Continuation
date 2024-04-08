function [psi2m,qm,wm,phim,psi2,q,w,phi] =  truncvec_2_specmat(U_t,Tq_vs_psi2,D2,D,Omega,rr,riprime,Nc,K)
len = (Nc+1)*(2*K+1);

psi2 = reshape(U_t(1:len),(Nc+1),(2*K+1));
phi    = reshape(U_t(len+1:2*len),(Nc+1),(2*K+1));


[Phi0,Phi0_r,W0] = base_state_flow(Omega, rr, riprime);
Phi0 = Phi0*ones(1,2*K+1);
%Utheta0 = -Phi0_r*ones(1,2*K+1);
W0 = W0*ones(1,2*K+1);
%Flux0 = -Phi0_r;
%
Phi = Phi0 + phi;

% to get initial data for charge density
psi2m = my_fft(psi2,2*K+1);
for m=0:1:K              %  Assign \hat{q}_m in accordance with the values of  \hat{psi2}_m.
    qm(:, m+1)=Tq_vs_psi2(:,:, m+1)*psi2m(:, m+1);
end
q = my_ifft(qm,2*K+1);


Phim = my_fft(Phi,2*K+1);
phim = my_fft(phi,2*K+1);

for m=0:1:K
    % this is the Laplacian
    Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2));
   wm(:,m+1) = -Lop*phim(:,m+1);
end
w = my_ifft(wm,2*K+1);



end


