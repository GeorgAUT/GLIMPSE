function [out]=cubic_nls_symm_Lawson(u_hat,mu,h)
%% Implementation of the symmetric Lawson method from
% Celledoni et al. Found Comput Math (2008) Prop. 2.2
%
% Input:    u_hat...(\hat{u}_n)_{n=-N/2+1}^{N/2}
%           mu...parameter in NLS
%           tau...stepsize

N=max(size(u_hat));
zeromode=zeros(size(u_hat));
zeromode(N/2)=1;

Nonlinearity=@(v) fftpi(abs(ifftpi(v)).^2.*ifftpi(v));

f = @(v) expilaplacian(u_hat,h)-i*mu*h*expilaplacian(Nonlinearity((expilaplacian(u_hat,h/2)+expilaplacian(v,-h/2))/2),h/2);

% Timestep using fixed-point iterations

solution=u_hat;
updatesize=1;
j=0;
while updatesize>1e-13 && j<100
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution);
end

out=solution;

end