function [output] = SME_backward_differ_semi_implicit_second_order(Tn, Tnp1, h,tol,maxit)
%% Second order semi-implicit scheme introduced by Xie et al. J. Comput. Phys. (2020)
% Scheme A
%
%  Input:   Tn.....value of T at time tn
%           Tnp1...value of T at time tn+h
%           h......time step
%           tol....required tolerance in GMRES
%           maxit..max number of iterations in GMRES

N=max(size(Tn));

% Auxilliary values and function for semi-implicit equation
intermediate=2*Tnp1-Tn;
normintermediate=sqrt(pagedot(intermediate,intermediate));
unitintermediate=intermediate./normintermediate;

g=@(x) 3/2*x-h*pagecross(intermediate,ifftpi_mat(dx_mat(dx_mat(fftpi_mat(x)))));
RHS=2*Tnp1-1/2*Tn;

% Analytical preconditioning
gvec=@(x) reshape(g(reshape(x,[1,3,N])),[3*N,1]);
gprecond=@(x) 2/3*pagedot(unitintermediate,x).*unitintermediate+1/h./normintermediate.*ifftpi_mat(dx_inv_mat(dx_inv_mat(fftpi_mat(pagecross(unitintermediate,x)))));
gprecondvec=@(x) reshape(gprecond(reshape(x,[1,3,N])),[3*N,1]);

% Choose preconditioning in the correct regime and solve implicit equations
if h<1/N
    output=gmres(gvec,reshape(RHS,[3*N,1]),[],tol,maxit);
else
    output=gmres(@(x) gprecondvec(gvec(x)),gprecondvec(reshape(RHS,[3*N,1])),[],tol,maxit);
end

Tnp2=reshape(output,[1,3,N]);

% Normalise and update the timestep values

output=Tnp2./(pagemtimes(Tnp2,pagetranspose(Tnp2))).^(1/2);
end

