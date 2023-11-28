function [out]=nls_nonlinear_part_SME(u,h)
%% Compute exact forward flow of partial_xu=-i/2|u|^{2}u
%
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
N=max(size(u));
u_phys=ifftpi(u);
out=fftpi(exp(i*h/2.*abs(u_phys).^2).*u_phys);
end