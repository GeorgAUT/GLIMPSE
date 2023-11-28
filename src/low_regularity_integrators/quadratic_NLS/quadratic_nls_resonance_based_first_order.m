function [phi1_out]=quadratic_nls_resonance_based_first_order(u_hat,mu,h)
%% Compute resonance-based first order scheme for quadratic NLS
% From Ostermann & Schratz 2018
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep

N=max(size(u_hat));
zeromode=zeros(size(u_hat));
zeromode(N/2)=1;

phi1_out=(1-2*i*mu*h*u_hat(N/2)).*expilaplacian(u_hat,h)+i*mu*zeromode*h*u_hat(N/2)^2+...
    mu/2*conv1(expilaplacian(dx_inv(u_hat),h),expilaplacian(dx_inv(u_hat),h))-...
    mu/2*expilaplacian(conv1(dx_inv(u_hat),dx_inv(u_hat)),h);
end