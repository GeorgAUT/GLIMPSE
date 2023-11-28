function [phi1_out]=kdv_resonance_based_first_order(u_hat,h)
%% Scheme (10) from Hofmanova & Schratz Numer. Math. (2017)
% 
% Input:    u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           h...timestep
phi1_out=expmdx3(u_hat,h)+1/6*conv1(expmdx3(dx_inv(u_hat),h),expmdx3(dx_inv(u_hat),h))...
    -1/6*expmdx3(conv1(dx_inv(u_hat),dx_inv(u_hat)),h);
end