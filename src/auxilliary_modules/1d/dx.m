function [dxu]=dx(v)
%% Compute spatial derivative
%
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}
N=max(size(v));
dxu=i*(-N/2+1:N/2)'.*v;
end
