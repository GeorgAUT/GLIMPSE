function [solution]=quadratic_nls_resonance_based_first_order_adjoint(u,mu,h)
%% Compute adjoint of resonance-based first order scheme for quadratic NLS
% From Ostermann & Schratz 2018
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u));
zeromode=zeros(size(u));
zeromode(N/2)=1;


% Compute adjoint step using fixed point iterations
f = @(v) expilaplacian(u,h)-2*i*mu*h*v(N/2).*v+i*mu*zeromode*h*v(N/2)^2-...
            mu/2*expilaplacian(conv1(expilaplacian(dx_inv(v),-h),expilaplacian(dx_inv(v),-h)),h)+...
            mu/2*conv1(dx_inv(v),dx_inv(v));       

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