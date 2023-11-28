function [dxu]=dx2d(v)
%% Compute spatial derivative in \partial_x direction
%
% Input:    u...(\hat{u}_{m,n)_{m,n=-N/2+1}^{N/2}
N=max(size(v));
dxu=1i*(-N/2+1:N/2).*v;
end
