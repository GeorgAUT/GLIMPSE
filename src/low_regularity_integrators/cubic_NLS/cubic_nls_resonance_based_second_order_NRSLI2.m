function [solution]=cubic_nls_resonance_based_second_order_NRSLI2(u_hat,mu,tau)
%% Non-resonant second-order low-regularity integrator
% based on arXiv:2302.00383 NRSLI2
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           tau...timestep

scale=1/2;

f = @(v) expilaplacian(u_hat,tau)-mu*tau/2*i*expilaplacian(...
   conv1(conv1(expilaplacian(u_hat,tau*scale),expilaplacian(u_hat,tau*scale)),expilaplacian(phi1laplacian(fftpi(conj(ifftpi(u_hat))),-2*tau),tau*scale))...
   +conv1(conv1(expilaplacian(v,-tau*scale),expilaplacian(v,-tau*scale)),expilaplacian(phi1laplacian(fftpi(conj(ifftpi(v))),-2*tau),tau*(scale+1)))...
   ,tau/2);

% Find update using fixed point iterations
solution=u_hat;
updatesize=1;
j=0;
while updatesize>tau^3 && j<150
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution); % Stopping criterion
end
end