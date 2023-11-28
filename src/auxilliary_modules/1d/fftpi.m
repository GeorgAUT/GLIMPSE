function [uhat1]=fftpi(u)
%% FFT on domain [-pi,pi]
% takes u=[u(2pi/N*(-N/2+1)),u(2pi/N*(-N/2+2)),...,u(2pi/N*(N/2))]
% maps to uhat=[\hat{u}_{-N/2+1},...,\hat{u}_{N/2}]
N=max(size(u));

u1=[u(N/2:end);u(1:N/2-1)]; % reordering the vector
uhat=fft(u1)/N; % Built-in FFT
uhat1=[uhat(N/2+2:end);uhat(1:N/2+1)]; % reordering
end