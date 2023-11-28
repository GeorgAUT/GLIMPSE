function [output] = SME_backward_differ_semi_implicit_first_order(Tn,h,tol,maxit)
%% First order semi-implicit scheme based on backward Euler in the spirit of
% E and Wang SIAM J. Numer. Anal. (2001)
%
%  Input:   Tn.....value of T at time tn
%           h......time step
%           tol....required tolerance in GMRES
%           maxit..max number of iterations in GMRES

N=max(size(Tn));

% Auxilliary values and function for semi-implicit equation
g=@(x) x-h*pagecross(Tn,ifftpi_mat(dx_mat(dx_mat(fftpi_mat(x)))));
RHS=Tn;


% Analytical preconditioning
gvec=@(x) reshape(g(reshape(x,[1,3,N])),[3*N,1]);
gprecond=@(x) pagedot(Tn,x).*Tn+1/h*ifftpi_mat(dx_inv_mat(dx_inv_mat(fftpi_mat(pagecross(Tn,x)))));
gprecondvec=@(x) reshape(gprecond(reshape(x,[1,3,N])),[3*N,1]);

% Choose preconditioning in the correct regime and solve implicit equations
if h<1/N
    output=gmres(gvec,reshape(RHS,[3*N,1]),[],tol,maxit);
else
    output=gmres(@(x) gprecondvec(gvec(x)),gprecondvec(reshape(RHS,[3*N,1])),[],tol,maxit);
end

Tnp1=reshape(output,[1,3,N]);

% Normalise and update the timestep values

output=Tnp1./(pagemtimes(Tnp1,pagetranspose(Tnp1))).^(1/2);
end