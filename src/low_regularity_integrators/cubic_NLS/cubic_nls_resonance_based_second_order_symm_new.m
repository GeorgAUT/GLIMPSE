function [solution]=cubic_nls_resonance_based_second_order_symm_new(u,mu,tau)
%% Symmetric resonance-based scheme
% based on arXiv:2305.16737 (4.18)
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u));

% Compute one time step in adjoint method by fixed pt iteration
f = @(v) expilaplacian(u,tau)-i*mu*tau/8*expilaplacian(conv1(conv1(u+expilaplacian(v,-tau),u+expilaplacian(v,-tau)),phi1laplacian(fftpi(conj(ifftpi(u+expilaplacian(v,-tau)))),-2*tau)-phi2laplacian(fftpi(conj(ifftpi(u+expilaplacian(v,-tau)))),-2*tau)),tau)...
                             -i*mu*tau/8*conv1(conv1(expilaplacian(u,tau)+v,expilaplacian(u,tau)+v),phi2laplacian(fftpi(conj(ifftpi(expilaplacian(u,-tau)+expilaplacian(v,-2*tau)))),-2*tau));

solution=u;
updatesize=1;
j=0;
while updatesize>1e-12 && j<150 % Iterations
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution); % Stopping criterion
end
    
end