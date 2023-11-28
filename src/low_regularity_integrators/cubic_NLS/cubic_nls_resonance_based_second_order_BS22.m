function [out]=cubic_nls_resonance_based_second_order_BS22(un,h)
%% Compute Second order low-regularity integrator
% From Bruned & Schratz 2022 Forum Pi
%
% Input:    un...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           h...timestep

un_phys=ifftpi(un); % Values in physical space
un_conj=fftpi(conj(un_phys)); % Fourier coefficients of conjugate
aux1=expilaplacian(un,h); % Auxilliary value

out=expilaplacian(un,h)...
            -i*h*expilaplacian(fftpi(un_phys.^2.*ifftpi(phi1laplacian(un_conj,-2*h)-phi2laplacian(un_conj,-2*h))),h)...
            -i*h*fftpi(ifftpi(aux1).*ifftpi(aux1).*ifftpi(phi2laplacian(expilaplacian(un_conj,h),-2*h)))...
            -h^2/2*expilaplacian(fftpi(abs(un_phys).^4.*un_phys),h);

end