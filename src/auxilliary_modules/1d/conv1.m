function [ustarv]=conv1(u,v)
%% Compute convolution of two functions by FFT and pointwise multiplication
%
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}
%           v...(\hat{v}_n)_{n=-N/2+1}^{N/2}
ustarv=fftpi(ifftpi(u).*ifftpi(v));
end