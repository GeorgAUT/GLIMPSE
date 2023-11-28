function [phi2_out]=phi2laplacian(u,h)
%% Compute \varphi_2(h*laplacian)u
%
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}
%           h...timestep
N=max(size(u));
laplacian=-transpose((-N/2+1:N/2).^2);
phi1=(exp(i*laplacian*h)-1)./(i*laplacian*h); % Precompute varphi1

phi2_out=(exp(i*laplacian*h)-phi1)./(i*laplacian*h).*u;
phi2_out(N/2)=1/2*u(N/2); % Treat special case of zeroth Fourier mode
end