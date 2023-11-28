%% In this example we evaluate the performance of our methods as a function of epsilon
% The results can be plotted using plotting_fn_epsilon.m
%
clear all
rng(10) % Fix random see for reproducability
%% Parameters

N_test=2^8; % Number of Fourier modes
N_ref=2^9;

N=N_ref;
Xvec=pi*(-N/2+1:N/2)'/N*2;
indexvec=(-N/2+1:N/2)';
indexvec(N/2)=1;

%% Initial conditions of O(1)
% Low-regularity initial conditions
theta=1.5; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;

% Smooth initial conditions
% u0=@(x) cos(x)./(2+sin(x));
% u0_eval=u0(Xvec);
% u0_hat=fftpi(u0_eval);

% Equation we consider is i\partial_t u= -\Delta u +\epsilon u^2
%% tau=0.25
Jmax=5;
error_sym_res=zeros(1,Jmax);
error_res=zeros(1,Jmax);
error_lie=zeros(1,Jmax);
error_strang=zeros(1,Jmax);
epsilon_j=zeros(1,Jmax); % For collecting timestep sizes

s=1 % Sobolev norm of error
    for j=1:Jmax
    
    %% Compute time-stepping for decreasing epsilon
    Jmax-j
    epsilon=2^(-j+1);
    T=1/epsilon; %Final time
    
    %% Compute reference solution with small time-step & resonance-based method
    tau=1e-4;
    M=floor(T/tau);
    tau=T/M;
    u_hat=u0_hat; % Start time-integration from initial condition
    for m=1:M
        M-m
        u_hat=quadratic_nls_resonance_based_first_order_adjoint(quadratic_nls_resonance_based_first_order(u_hat,epsilon,tau/2),epsilon,tau/2);
    end
    
    u_hat_ref=u_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Save reference solution
    
    tau=0.25
    M=floor(T/tau);
    tau=T/M;
    u_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with standard resonance-based
    v_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with symmetric resonance-based
    w_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with lie splitting
    z_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with Strang splitting
    tn=0;
    for m=1:M
        % Resonance-based methods
        u_hat=quadratic_nls_resonance_based_symm_second_order(u_hat,epsilon,tau);
    
        v_hat=quadratic_nls_resonance_based_first_order_adjoint(quadratic_nls_resonance_based_first_order(v_hat,epsilon,tau/2),epsilon,tau/2);
        % Strang splitting
    
        z_hat=expilaplacian(z_hat,tau/2);
        z_hat=nls_quadratic_nonlinear_part(z_hat,epsilon,tn,tau);
        z_hat=expilaplacian(z_hat,tau/2);
    
        % lie splitting
    
        w_hat=nls_quadratic_nonlinear_part(w_hat,epsilon,tn,tau);
        w_hat=expilaplacian(w_hat,tau);
        tn=tn+tau;
    end
    
    
    N=N_test;
    Xvec=pi*(-N/2+1:N/2)'/N*2;
    indexvec=(-N/2+1:N/2)';
    indexvec(N/2)=1;
    
    error_sym_res(j)=norm((u_hat-u_hat_ref).*abs(indexvec).^s);
    error_res(j)=norm((v_hat-u_hat_ref).*abs(indexvec).^s);
    error_lie(j)=norm((w_hat-u_hat_ref).*abs(indexvec).^s);
    error_strang(j)=norm((z_hat-u_hat_ref).*abs(indexvec).^s);
    epsilon_j(j)=epsilon;
    
    save(strcat('data/long_time_function_of_epsilon_tau_0_25_H1_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'_v2.mat'))
end

save(strcat('data/long_time_function_of_epsilon_tau_0_25_H1_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'_v2.mat'))