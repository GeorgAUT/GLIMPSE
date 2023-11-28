function [solution,Gamma_hat]=cubic_nls_crank_nicolson_relaxed(u,Gamma_hat,mu,tau)
%% Energy preserving relaxed Crank-Nicolson for cubic NLS
% as introduced by Besse 2004, SIAM J. Numer. Anal.
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u));

% Then take step in Gamma

Gamma_hat=-Gamma_hat+2*fftpi(abs(ifftpi(u)).^2);

% Define implicit system

g=@(v) v-i*tau/2*dx(dx(v))+i*tau/2*mu*conv1(Gamma_hat,v);

solution=gmres(g,u+i*tau/2*dx(dx(u))-i*mu*tau/2*conv1(Gamma_hat,u),100,1e-15,1e3);

    
end