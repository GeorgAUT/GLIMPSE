function [phi1_out]=kdv_resonance_based_second_order(u_hat,h)
%% Scheme (17) from Hofmanova & Schratz Numer. Math. (2017)
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           h...timestep
phi1_out=expmdx3(u_hat,h)+1/6*conv1(expmdx3(dx_inv(u_hat),h),expmdx3(dx_inv(u_hat),h))...
    -1/6*expmdx3(conv1(dx_inv(u_hat),dx_inv(u_hat)),h)...
    -1/18*dx_inv(conv1(expmdx3(dx_inv(dx_inv(u_hat)),h),expmdx3(dx_inv(conv1(u_hat,u_hat)),h)))...
    +1/18*expmdx3(dx_inv(dx_inv(conv1(dx_inv(dx_inv(u_hat)),dx_inv(conv1(u_hat,u_hat))))),h);

N=max(size(u_hat));
phi1_out(N/2)=0;
end