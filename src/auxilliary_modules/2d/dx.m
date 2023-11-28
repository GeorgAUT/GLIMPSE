function [dxu]=dx(v)
%% Compute spatial derivative in x direction
%
% Input:    u...(\hat{u}_{m,n)_{m,n=-N/2+1}^{N/2}
N=max(size(v));
dxu=1i*(-N/2+1:N/2).*v;
end
