function [out]=cubic_nls_2d_resonance_based_first_order(un,mu,h)
%% Compute resonance-based first order scheme in 2d
% From Ostermann & Schratz 2018
% 
% Input:    un...value of (\hat{u}_{m,n})_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
out=expilaplacian2d(un-mu*h*i*conv2d(un,conv2d(un,phi1laplacian2d(fft2pi(conj(ifft2pi(un))),-2*h))),h);
end