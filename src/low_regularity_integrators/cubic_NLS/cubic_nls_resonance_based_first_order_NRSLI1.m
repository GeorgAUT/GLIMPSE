function [out]=cubic_nls_resonance_based_first_order_NRSLI1(u_hat,mu,tau)
%% Non-resonant first-order low-regularity integrator
% based on arXiv:2302.00383 NRSLI1
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           tau...timestep

N=max(size(u_hat));
aux1=conv1(fftpi(conj(ifftpi(u_hat)))-phi1laplacian(fftpi(conj(ifftpi(u_hat))),-2*tau),u_hat);
laplacian=-transpose((-N/2+1:N/2).^2);
phi1=(exp(-i*2*laplacian*tau)-1)./(-i*2*laplacian*tau);
phi1(N/2)=1;
aux2=abs(u_hat).^2.*u_hat.*(1-phi1);

out=cubic_nls_resonance_based_first_order(u_hat,mu,tau)-i*mu*tau*expilaplacian(2*u_hat*aux1(N/2)-aux2,tau);
end