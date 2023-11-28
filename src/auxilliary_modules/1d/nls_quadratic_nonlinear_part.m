function [u_hat_out]=nls_quadratic_nonlinear_part(u_hat,epsilon,tn,tau)
%% Compute exact forward flow of partial_tu=-i*u^2
%
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
c=(1-i*tn*epsilon*ifftpi(u_hat))./ifftpi(u_hat);
c(find(abs(c)<1e-15))=1;

u_hat_out=fftpi(1./(c+i*epsilon*(tn+tau)));
end