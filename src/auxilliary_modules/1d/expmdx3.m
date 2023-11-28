function [phi1_out]=expmdx3(u,h)
%% Compute exact forward flow of partial_tu=-(partial_x)^{3}u
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}
N=max(size(u));
phi1_out=exp(i*h.*((-N/2+1:N/2).^3)').*u;
end