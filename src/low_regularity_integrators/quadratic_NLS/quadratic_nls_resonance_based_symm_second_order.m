function [solution]=quadratic_nls_resonance_based_symm_second_order(u,mu,h)
%% Symmetric second order low-regularity integrator
% based on arXiv:2302.00383 SLI2
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep

N=max(size(u));
zeromode=zeros(size(u));
zeromode(N/2)=1;

% Update is computed with fixed point iterations

f = @(v) expilaplacian(u,h)-i*mu*h*v(N/2).*v-i*mu*h*u(N/2).*expilaplacian(u,h)+...
            i*mu*zeromode*h/2*u(N/2)^2+i*mu*zeromode*h/2*v(N/2)^2-...
            mu/4*expilaplacian(conv1(expilaplacian(dx_inv(v),-h),expilaplacian(dx_inv(v),-h)),h)-...
            mu/4*expilaplacian(conv1(dx_inv(u),dx_inv(u)),h)+...
            mu/4*conv1(dx_inv(v),dx_inv(v))+...
            mu/4*conv1(expilaplacian(dx_inv(u),h),expilaplacian(dx_inv(u),h));

solution=u;
updatesize=1;
j=0;
while updatesize>h^3 && j<150
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution); % Stopping criterion
end
end