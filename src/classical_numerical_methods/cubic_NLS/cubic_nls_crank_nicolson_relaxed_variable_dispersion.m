function [solution,Gamma_hat]=cubic_nls_crank_nicolson_relaxed_variable_dispersion(u,Gamma_hat,mu,tau,dB)
%% A relaxed CN method for cubic NLSE with white noise dispersion based on
% Besse 2004, SIAM J. Numer. Anal.
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           tau...timestep
%           dB=B(t_{n+1})-B(t_n)
N=max(size(u));

% Then take step in Gamma

Gamma_hat=-Gamma_hat+2*fftpi(abs(ifftpi(u)).^2);

% Define implicit system

g=@(v) v-i*dB/2*dx(dx(v))+i*tau/2*mu*conv1(Gamma_hat,v);

[solution,flag]=gmres(g,u+i*dB/2*dx(dx(u))-i*mu*tau/2*conv1(Gamma_hat,u),12,1e-8,1e4);

    
end