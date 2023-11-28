function [solution]=kdv_onestageLawson1(u_hat,h)
%% Implementation of the symmetric Lawson method from
% Celledoni et al. Found Comput Math (2008) Prop. 2.2
%
% Input:    u_hat...(\hat{u}_n)_{n=-N/2+1}^{N/2}
%           h...stepsize


N=@(w) 1/2*dx(conv1(w,w));
f=@(v) N(expmdx3(u_hat,h/2)+h/2*v);

% Solve using fixed point iterations
% initial_guess=u_hat;
N1=u_hat;
j=0;
N1err=1;
while j<10 && N1err>max(1e-16,1e-2*h^3)
    N1old=N1;
    N1=f(N1);
    N1err=norm(N1-N1old)/norm(N1); % Stopping criterion
    j=j+1;
end
solution=expmdx3(u_hat,h)+h*expmdx3(N1,h/2);

end