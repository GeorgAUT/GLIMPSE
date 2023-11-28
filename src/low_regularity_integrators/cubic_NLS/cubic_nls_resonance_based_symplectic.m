function [out]=cubic_nls_resonance_based_symplectic(u_hat,mu,h)
%% Symplectic resonance-based scheme
% based on arXiv:2205.05024 Example 4.6
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
N=max(size(u_hat));
zeromode=zeros(size(u_hat));
zeromode(N/2)=1;

I1=@(v) i/2*dx_inv(expilaplacian(fftpi(ifftpi(expilaplacian(fftpi(ifftpi(v).^2),h))...
            .*ifftpi(dx_inv(expilaplacian(fftpi(conj(ifftpi(v))),-h)))),-h)...
            -fftpi(ifftpi(v).^2.*ifftpi(dx_inv(fftpi(conj(ifftpi(v)))))))...
            +h*zeromode.*fftpi(ifftpi(v).^2.*conj(ifftpi(v)))...
            -h*conj(v(N/2))*fftpi(ifftpi(v).^2).*zeromode+h*conj(v(N/2))*fftpi(ifftpi(v).^2);
    
I2=@(v) i/2*fftpi(conj(ifftpi(v)).*(ifftpi(expilaplacian(fftpi(ifftpi(dx_inv(expilaplacian(v,h))).^2),-h))...
                -ifftpi(dx_inv(v)).^2))...
                +2*h*v(N/2)*fftpi(abs(ifftpi(v)).^2)-h*v(N/2)^2*fftpi(conj(ifftpi(v)));
            
            
I3=@(v) (-h)*fftpi(ifftpi(v).^2.*conj(ifftpi(v)));


% Construct function for fixed point iterations

midpt=@(v) (u_hat+expilaplacian(v,-h))/2;

f = @(v) expilaplacian(u_hat,h)-i*mu*expilaplacian(I1(midpt(v))+I2(midpt(v))+I3(midpt(v)),h);

% Fixed point iterations

solution=u_hat;
updatesize=1;
j=0;
while updatesize>1e-16 && j<400 %h^4 && j<100
    solution_old=solution;
    solution=f(solution_old);
    j=j+1;
    updatesize=norm(solution_old-solution);
end

out=solution;

end