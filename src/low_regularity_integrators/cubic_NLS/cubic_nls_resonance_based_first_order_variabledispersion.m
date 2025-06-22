function [out]=cubic_nls_resonance_based_first_order_variabledispersion(un,mu,psin,taulist)
%% Resonance-based first order scheme for NLSE with variable dispersion
% From Cui & Maierhofer 2025 (arXiv:2503.19346)
% 
% Input:    un...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           mu...parameter in nonlinearity of NLS
%           psin...Brownian motion increment
%           taulist...list of subintervals for Duhamel discretisation

% Some parameters specific to Wong-Zakai
psintau=psin(end)-psin(1);
taulist=taulist(2:end)-taulist(1:end-1);
c=(psin(2:end)-psin(1:end-1))./taulist; %slope of linear interpolant locally
b=psin(2:end)-c.*taulist; % Shift of interpolant locally


% Compute terms as in the deterministic case
out1=expilaplacian(un,psintau);
unconj=fftpi(conj(ifftpi(un)));

% Carefully include the linear interpolants in the exponent
aux=sum(taulist.*expilaplacian(phi1laplacian(unconj,-2*c.*taulist),-2*b),2);

out2=-mu*i*conv1(un,conv1(un,aux));
out2=expilaplacian(out2,psintau);

out=out1+out2;

end