%% Comparison of convergence behaviour of new methods
% To plot the results please use plotting_fn_tau.m
% To compute reference solution, use compute_reference_solution.m
%
clear all
rng(10) % Fix random see for reproducability
%% Parameters

N=2^8; % Number of Fourier modes
Xvec=pi*(-N/2+1:N/2)'/N*2;
indexvec=(-N/2+1:N/2)';
indexvec(N/2)=1;

%% Initial conditions of O(1)
% Low-regularity initial conditions
theta=2; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;
% u0_hat(N/2)=0;
% Smooth initial conditions
% u0=@(x) cos(x)./(2+sin(x));
% u0_eval=u0(Xvec);
% u0_hat=fftpi(u0_eval);

%% Equation we consider is i\partial_t u= -\Delta u +\epsilon^2 |u|^2u
% Epsilon 0.1

epsilon=0.5;
T=1/epsilon^2; %Final time

%% Compute reference solution with small time-step & resonance-based method
load('data/reference_solution_epsilon_0-5_H2_N_256.mat','u0_hat','u_hat_ref','N_ref')
u_hat_ref=u_hat_ref(N_ref/2-N/2+1:N_ref/2+N/2); % Save reference solution
u0_hat=u0_hat(N_ref/2-N/2+1:N_ref/2+N/2);

%% Compute time-stepping for decreasing tau
Jmax=20;
error_sym_res=zeros(1,Jmax);
error_res=zeros(1,Jmax);
error_strang=zeros(1,Jmax);
error_lie=zeros(1,Jmax);
error_res_long_time=zeros(1,Jmax);
error_res_long_time_second_order=zeros(1,Jmax);
error_res_long_time_second_order_yue=zeros(1,Jmax);
tau_j=zeros(1,Jmax); % For collecting timestep sizes

s=1 % Sobolev norm of error
for ll=0:Jmax/8-1
    parfor j=ll*8+1:ll*8+8
    Jmax-j
    %M=400*j+1000;
    tau=1/(5*j)^1.0;
    M=floor(T/tau);
    tau=T/M;

    u_hat=u0_hat; % Sln computed with symmetric resonance-based
    v_hat=u0_hat; % Sln computed with standard resonance-based
    w_hat=u0_hat; % Sln computed with strang splitting
    zz_hat=u0_hat; % Sln computed with lie splitting
    zzz_hat=u0_hat; % Sln computed with new resonance_longtime
    zzzz_hat=u0_hat; % Sln computed with new resonance_longtime second order
    zzzzz_hat=u0_hat; % Sln computed with new resonance_longtime second order

    for m=1:M
        %M-m
        % Resonance-based methods
        
        u_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(u_hat,epsilon^2,tau/2),epsilon^2,tau/2);
        v_hat=cubic_nls_resonance_based_first_order(v_hat,epsilon^2,tau);
        
        zzz_hat=cubic_nls_resonance_based_first_order_NRSLI1(zzz_hat,epsilon^2,tau);

        zzzzz_hat=cubic_nls_resonance_based_second_order_NRSLI2(zzzzz_hat,epsilon^2,tau);%
        zzzz_hat=cubic_nls_resonance_based_first_order_NRSLI1_adjoint(cubic_nls_resonance_based_first_order_NRSLI1(zzzz_hat,epsilon^2,tau/2),epsilon^2,tau/2);
        
        % Strang splitting

        w_hat=expilaplacian(w_hat,tau/2);
        w=ifftpi(w_hat);
        w_hat=fftpi(exp(-i*tau*epsilon^2*abs(w).^2).*w);
        w_hat=expilaplacian(w_hat,tau/2);
        
        
        % Lie splitting

        zz_hat=expilaplacian(zz_hat,tau);
        zz=ifftpi(zz_hat);
        zz_hat=fftpi(exp(-i*tau*epsilon^2*abs(zz).^2).*zz);
    end

    error_sym_res(j)=norm((u_hat-u_hat_ref).*abs(indexvec).^s);
    error_res(j)=norm((v_hat-u_hat_ref).*abs(indexvec).^s);
    error_strang(j)=norm((w_hat-u_hat_ref).*abs(indexvec).^s);
    error_lie(j)=norm((zz_hat-u_hat_ref).*abs(indexvec).^s);
    error_res_long_time(j)=norm((zzz_hat-u_hat_ref).*abs(indexvec).^s);
    error_res_long_time_second_order(j)=norm((zzzz_hat-u_hat_ref).*abs(indexvec).^s);
    error_res_long_time_second_order_yue(j)=norm((zzzzz_hat-u_hat_ref).*abs(indexvec).^s);
    tau_j(j)=tau;
    end

    save(strcat('data/long_time_function_of_tau_epsilon_0-5_H',num2str(theta),'_data_N_',num2str(N),'_v2.mat'))

end

save(strcat('data/long_time_function_of_tau_epsilon_0-5_H',num2str(theta),'_data_N_',num2str(N),'_v2.mat'))
