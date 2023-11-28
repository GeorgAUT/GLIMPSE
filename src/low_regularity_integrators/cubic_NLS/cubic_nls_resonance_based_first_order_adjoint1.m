function [solution]=cubic_nls_resonance_based_first_order_adjoint1(u,mu,tau)
%% Compute adjoint of resonance-based first order scheme
% based on Ostermann & Schratz 2018
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u));


% Compute one time step in adjoint method by fixed pt iteration
f = @(v) expilaplacian(u,tau)-mu*tau*i*conv1(v,conv1(v,phi1laplacian(fftpi(conj(ifftpi(v))),2*tau))); % Function for fixed point iterations

solution=u;
updatesize=1;
j=0;
while updatesize>1e-14 && j<150 % Iterations
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution); % Stopping criterion
end
    
end