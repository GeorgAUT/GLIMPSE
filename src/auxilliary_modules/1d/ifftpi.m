function [u1]=ifftpi(uhat)
%% IFFT on domain [-pi,pi]
% takes uhat=[\hat{u}_{-N/2+1},...,\hat{u}_{N/2}]
% maps to u=[u(2pi/N*(-N/2+1)),u(2pi/N*(-N/2+2)),...,u(2pi/N*(N/2))]
N=max(size(uhat));


uhat1=[uhat(N/2:end);uhat(1:N/2-1)]; % reordering the vector
u=ifft(uhat1)*N; % Built-in IFFT and rescaling
u1=[u(N/2+2:end);u(1:N/2+1)]; % reordering
end