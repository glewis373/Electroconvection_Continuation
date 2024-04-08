function [psi2m_new,qm_new,Phim_new,Wm_new,flux_new] = second_order_time_step(psi2m_new,...
    qm_old,qm_now,qm_new,Phim_new,Wm_old,Wm_now,Wm_new,flux_now,flux_old,flux_nonl_now,...
    flux_nonl_old,Jqphi_hat_old,Jqphi_hat_now,...
    Jwphi_hat_old,Jwphi_hat_now,Jpsi2q_hat_old,Jpsi2q_hat_now,D,D2,K,riprime,...
    dt,Tq_vs_psi2,Nc,T_inv,c_int,wm1,wm3,phim1,phim3,Mi,phi_BC,Ra,Pr,in_speed,It)


flux_new = second_order_flux_update(flux_now,flux_old,flux_nonl_now,...
        flux_nonl_old,riprime,in_speed,Pr,It,D,D2,dt,Nc);

for m=0:1:K % Loop for each Fourier mode
    % Lop_PT = 4*D2 + 2*riprime_invM.*D - m^2* full(sparse(1:Nc+1, 1:Nc+1, 1./riprime_square )); % Laplace operator matrix (checked)
    Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2));
    
    % first step: solve for the new surface charge and the new 2d
    % electric potential
    [psi2m_new,qm_new] = second_order_charge_update(psi2m_new,qm_new,qm_now,qm_old,dt,...
            Tq_vs_psi2,Lop,m,Nc,Jqphi_hat_now,Jqphi_hat_old);

    % second step: now solve for the new vorticity and fluid potential
        [Wm_new,Phim_new] = second_order_vorticity_update(Wm_new,Phim_new,Wm_now,Wm_old,flux_new,D,...
	  wm1,wm3,phim1,phim3,Mi,phi_BC,Ra,Pr,in_speed,dt,It,Lop,m,Nc,T_inv,c_int,K,...
            Jwphi_hat_now,Jwphi_hat_old,Jpsi2q_hat_now,Jpsi2q_hat_old);
end
        