function [solution]=cubic_nls_crank_nicolson(u,mu,tau)
%% Classical Crank-Nicolson scheme
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u));

% Compute one time step by fixed pt iteration
g=@(v1,v0) i*dx(dx(v1+v0))/2-i*mu/4*fftpi((abs(ifftpi(v1)).^2+(abs(ifftpi(v0)).^2)).*ifftpi(v1+v0));

f = @(v) u+tau*(g(u,v));

solution=u;
updatesize=1;
j=0;
while updatesize>1e-14 && j<450 % Iterations
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution); % Stopping criterion
    if isnan(sum(solution))
        break;
    end
end
    
end