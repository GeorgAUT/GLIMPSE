function [uhat]=fft2pi(u)
%% FFT on domain [-pi,pi]x[-pi,pi]
% takes u=u(Xvec,Yvec)
% maps to uhat=[\hat{u}_{-N/2+1,-N/2+1},...,\hat{u}_{N/2,-N/2+1};
%               ...;
%               \hat{u}_{-N/2+1,N/2},...,\hat{u}_{N/2,N/2}]
N=max(size(u));

u=[u(:,N/2:end),u(:,1:N/2-1)]; % reordering the vector
u=[u(N/2:end,:);u(1:N/2-1,:)]; % reordering the vector
uhat=fft2(u)/N/N; % Built-in FFT
uhat=[uhat(:,N/2+2:end),uhat(:,1:N/2+1)]; % reordering
uhat=[uhat(N/2+2:end,:);uhat(1:N/2+1,:)]; % reordering
end