function [solution]=cubic_nls_2d_resonance_based_second_order_symm_new(u,mu,tau)
%% Symmetric resonance-based scheme
% based on arXiv:2305.16737 (4.18)
% 
% Input:    u...value of (\hat{u}_{m,n})_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u));

% Compute one time step in method by fixed pt iteration
f = @(v) expilaplacian2d(u,tau)-i*mu*tau/8*expilaplacian2d(conv2d(conv2d(u+expilaplacian2d(v,-tau),u+expilaplacian2d(v,-tau)),phi1laplacian2d(fft2pi(conj(ifft2pi(u+expilaplacian2d(v,-tau)))),-2*tau)-phi2laplacian2d(fft2pi(conj(ifft2pi(u+expilaplacian2d(v,-tau)))),-2*tau)),tau)...
                             -i*mu*tau/8*conv2d(conv2d(expilaplacian2d(u,tau)+v,expilaplacian2d(u,tau)+v),phi2laplacian2d(fft2pi(conj(ifft2pi(expilaplacian2d(u,-tau)+expilaplacian2d(v,-2*tau)))),-2*tau)); % Function for fixed point iterations

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