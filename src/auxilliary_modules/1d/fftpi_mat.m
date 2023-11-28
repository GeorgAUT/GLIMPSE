function [uhat1]=fftpi_mat(u)
%% FFT on domain [-pi,pi]
% takes u=[u(2pi/N*(-N/2+1)),u(2pi/N*(-N/2+2)),...,u(2pi/N*(N/2))]
% maps to uhat=[\hat{u}_{-N/2+1},...,\hat{u}_{N/2}]
% Applies to vectors where the final coordinate corresponds to length
N=max(size(u,3));
% reordering the vector
%uvec=u0_cts(pi*(-N/2+1:N/2)/N*2);
u1=cat(3,u(:,:,N/2:end),u(:,:,1:N/2-1));
uhat=permute(fft(permute(u1,[3,1,2]))/N,[2,3,1]);
uhat1=cat(3,uhat(:,:,N/2+2:end),uhat(:,:,1:N/2+1));
end