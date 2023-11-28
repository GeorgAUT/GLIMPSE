function [out]=cubic_nls_resonance_based_first_order(un,mu,h)
%% Compute resonance-based first order scheme
% From Ostermann & Schratz 2018
% 
% Input:    un...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           h...timestep
out=expilaplacian(un-mu*h*i*conv1(un,conv1(un,phi1laplacian(fftpi(conj(ifftpi(un))),-2*h))),h);
end