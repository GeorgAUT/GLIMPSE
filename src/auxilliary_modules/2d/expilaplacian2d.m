function [out]=expilaplacian2d(u,h)
%% Compute exact forward flow of partial_tu=i(\Delta)^{2}u in 2D
% Input:    u=[\hat{u}_{-N/2+1,-N/2+1},...,\hat{u}_{N/2,-N/2+1};
%               ...;
%               \hat{u}_{-N/2+1,N/2},...,\hat{u}_{N/2,N/2}]

N=max(size(u));
out=exp(-i*h.*(((-N/2+1:N/2).^2)'+(-N/2+1:N/2).^2)).*u;
end