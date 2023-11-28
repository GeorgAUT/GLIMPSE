function [u]=ifft2pi(uhat)
%% IFFT on domain [-pi,pi]
% takes uhat=[\hat{u}_{-N/2+1,-N/2+1},...,\hat{u}_{N/2,-N/2+1};
%               ...;
%               \hat{u}_{-N/2+1,N/2},...,\hat{u}_{N/2,N/2}]
%   maps to     u=u(Xvec,Yvec)
% 
N=max(size(uhat));


uhat=[uhat(:,N/2:end),uhat(:,1:N/2-1)]; % reordering the vector
uhat=[uhat(N/2:end,:);uhat(1:N/2-1,:)]; % reordering the vector
u=ifft2(uhat)*N*N; % Built-in IFFT and rescaling
u=[u(:,N/2+2:end),u(:,1:N/2+1)]; % reordering
u=[u(N/2+2:end,:);u(1:N/2+1,:)]; % reordering
end