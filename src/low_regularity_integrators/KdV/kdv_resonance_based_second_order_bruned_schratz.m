function [phi1_out]=kdv_resonance_based_second_order_bruned_schratz(u_hat,h)
%% Scheme (151) from Bruned & Schratz 2022 Forum Pi
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           h...timestep

phi1_out=expmdx3(u_hat,h)+1/6*conv1(expmdx3(dx_inv(u_hat),h),expmdx3(dx_inv(u_hat),h))...
    -1/6*expmdx3(conv1(dx_inv(u_hat),dx_inv(u_hat)),h)...
    +h^2/4*expmdx3(phi1laplacian(dx(conv1(u_hat,dx(conv1(u_hat,u_hat)))),h),h);

N=max(size(u_hat));
phi1_out(N/2)=0;
end