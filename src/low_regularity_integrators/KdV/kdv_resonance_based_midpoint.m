function [solution]=kdv_resonance_based_midpoint(u_hat,h)
%% New resonance-based midpoint scheme from Example 4.6 in arXiv:2205.05024
% 
% % Input:  u_hat...value of (\hat{u}_n)_{n=-N/2+1}^{N/2} at time t_n
%           h...timestep
f = @(v) expmdx3(u_hat,h)+1/6*conv1(expmdx3(dx_inv((expmdx3(v,-h)+u_hat)/2),h),expmdx3(dx_inv((expmdx3(v,-h)+u_hat)/2),h))...
-1/6*expmdx3(conv1(dx_inv((expmdx3(v,-h)+u_hat)/2),dx_inv((expmdx3(v,-h)+u_hat)/2)),h);

% Time step is computed using fixedpoint iterations
solution=u_hat;
j=0;
err=1;
while j<150 & err>1e-16
    j=j+1;
    solution_old=solution;
    solution=f(solution);
    err=norm(solution-solution_old)/norm(solution);
end
end