function [out]=expilaplacian(u,h)
%% Compute exact forward flow of partial_tu=i(partial_x)^{2}u
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}

N=max(size(u));
out=exp(-i*h.*((-N/2+1:N/2).^2)').*u;
end