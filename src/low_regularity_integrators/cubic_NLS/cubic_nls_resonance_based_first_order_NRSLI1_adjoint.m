function [solution]=cubic_nls_resonance_based_first_order_NRSLI1_adjoint(u,mu,tau)
%% Adjoint of non-resonant first-order low-regularity integrator
% based on arXiv:2302.00383 NRSLI1
% 
% Input:    u...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           tau...timestep

N=max(size(u));
a=zeros(1,N); a(N/2)=1;
aux1=@(u,tau) a*(conv1(fftpi(conj(ifftpi(u)))-phi1laplacian(fftpi(conj(ifftpi(u))),-2*tau),u));
laplacian=-transpose((-N/2+1:N/2).^2);
phi1=(exp(-i*2*laplacian*(-tau))-1)./(-i*2*laplacian*(-tau));
phi1(N/2)=1;
aux2=@(u,tau) abs(u).^2.*u.*(1-phi1);

% Find value using fixed-point iterations
f = @(v) expilaplacian(u,tau)-mu*tau*i*conv1(v,conv1(v,phi1laplacian(fftpi(conj(ifftpi(v))),2*tau)))...
    -i*mu*tau*(2*v*aux1(v,-tau)-aux2(v,-tau));

solution=u;
updatesize=1;
j=0;
while updatesize>tau^4 && j<150
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution); % Stopping criterion
end

end