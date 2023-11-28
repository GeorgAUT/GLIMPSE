function [out]=nls_nonlinear_part2d(u,h)
%% Compute exact forward flow of partial_tu=-i|u|^{2}u in 2D
%
% Input:    u...value of (\hat{u}_{m,n})_{m,n=-N/2+1}^{N/2} at time t_n
N=max(size(u));
u_phys=ifft2pi(u);
out=fft2pi(exp(-i*h.*abs(u_phys).^2).*u_phys);
end