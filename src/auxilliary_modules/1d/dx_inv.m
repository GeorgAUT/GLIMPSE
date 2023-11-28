function [dxinvu]=dx_inv(v)
%% Compute inverse of spatial derivative
%
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}
N=max(size(v));
dxinvu=(i*(-N/2+1:N/2)').^(-1).*v;
dxinvu(N/2)=0; % Set \hat{u}_0=0
end
