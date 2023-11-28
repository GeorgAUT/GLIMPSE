function [out]=cubic_nls_2d_resonance_based_second_order_BS22(un,mu,h)
%% Compute Second order low-regularity integrator in 2D
% From Bruned & Schratz 2022
%
% Input:    un...value of (\hat{u}_{m,n})_{n=-N/2+1}^{N/2} at time t_n
%           h...timestep

un_phys=ifft2pi(un); % Values in physical space
un_conj=fft2pi(conj(un_phys)); % Fourier coefficients of conjugate
aux1=expilaplacian2d(un,h); % Auxilliary value

out=expilaplacian2d(un,h)...
            -i*mu*h*expilaplacian2d(fft2pi(un_phys.^2.*ifft2pi(phi1laplacian2d(un_conj,-2*h)-phi2laplacian2d(un_conj,-2*h))),h)...
            -i*mu*h*fft2pi(ifft2pi(aux1).*ifft2pi(aux1).*ifft2pi(phi2laplacian2d(expilaplacian2d(un_conj,h),-2*h)))...
            -h^2*mu^2/2*expilaplacian2d(fft2pi(abs(un_phys).^4.*un_phys),h);

end