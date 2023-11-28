function [solution]=cubic_nls_2d_resonance_based_first_order_adjoint(u,mu,tau)
%% Compute adjoint of resonance-based first order scheme in 2d
% based on Ostermann & Schratz 2018
% 
% Input:    u...value of (\hat{u}_{m,n})_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           tau...timestep
N=max(size(u));


% Compute one time step in adjoint method by fixed pt iteration
f = @(v) expilaplacian2d(u,tau)-mu*tau*1i*conv2d(v,conv2d(v,phi1laplacian2d(fft2pi(conj(ifft2pi(v))),2*tau))); % Function for fixed point iterations

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