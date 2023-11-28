function [dxu]=dx_mat(v)
%% Compute spatial derivative
%
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}
% Applies to vectorfields where the final coordinate corresponds to length
N=max(size(v,3));
Fourier=zeros(1,1,N);
Fourier(1,1,:)=i*(-N/2+1:N/2);
dxu=Fourier.*v;
end
