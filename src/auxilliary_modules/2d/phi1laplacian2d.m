function [phi1_out]=phi1laplacian2d(u,h)
%% Compute \varphi_1(h*laplacian)u in 2D
%
% Input:    u...(\hat{u}_{m,n})_{m,n=-N/2+1}^{N/2}
%           h...timestep
N=max(size(u));
laplacian=-transpose((-N/2+1:N/2).^2)-(-N/2+1:N/2).^2;
phi1_out=(exp(i*laplacian*h)-1)./(i*laplacian*h).*u;
phi1_out(N/2,N/2)=u(N/2,N/2); % Note the special case for zeroth mode
end