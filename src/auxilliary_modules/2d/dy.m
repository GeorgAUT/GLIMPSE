function [dyu]=dy(v)
%% Compute spatial derivative in y direction
%
% Input:    u...(\hat{u}_{m,n)_{m,n=-N/2+1}^{N/2}
N=max(size(v));
dyu=1i*(-N/2+1:N/2)'.*v;
end
