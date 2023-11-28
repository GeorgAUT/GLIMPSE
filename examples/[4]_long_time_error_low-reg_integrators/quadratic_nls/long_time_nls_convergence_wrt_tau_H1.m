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
theta=1.5; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;

% Smooth initial conditions
% u0=@(x) cos(x)./(2+sin(x));
% u0_eval=u0(Xvec);
% u0_hat=fftpi(u0_eval);

% Equation we consider is i\partial_t u= -\Delta u +\epsilon u^2
%% Epsilon=0.1
epsilon=0.1;
T=1/epsilon; %Final time

%% Compute reference solution with small time-step & resonance-based method


load('data/reference_solution_epsilon_0-1_H1_N_256.mat','u0_hat','u_hat_ref','N_ref')
u_hat_ref=u_hat_ref(N_ref/2-N/2+1:N_ref/2+N/2); % Save reference solution
u0_hat=u0_hat(N_ref/2-N/2+1:N_ref/2+N/2);
%% Compute time-stepping for decreasing tau
Jmax=20;
error_sym_res=zeros(1,Jmax);
error_res=zeros(1,Jmax);
error_lie=zeros(1,Jmax);
error_strang=zeros(1,Jmax);
tau_j=zeros(1,Jmax); % For collecting timestep sizes

s=1 % Sobolev norm of error
for ll=0:Jmax/10-1
    parfor j=ll*10+1:ll*10+10
    Jmax-j
    tau=1/(5*j);
    M=floor(T/tau);
    tau=T/M;
    %M=1
    u_hat=u0_hat; % Sln computed with symmetric resonance-based
    v_hat=u0_hat; % Sln computed with standard resonance-based
    w_hat=u0_hat; % Sln computed with lie splitting
    z_hat=u0_hat; % Sln computed with lie splitting
    tn=0;
    for m=1:M
        % Resonance-based methods
        u_hat=quadratic_nls_resonance_based_symm_second_order(u_hat,epsilon,tau);

        v_hat=quadratic_nls_resonance_based_first_order(v_hat,epsilon,tau);
        % Strang splitting
        
        z_hat=expilaplacian(z_hat,tau/2);
        z_hat=nls_quadratic_nonlinear_part(z_hat,epsilon,tn,tau);
        z_hat=expilaplacian(z_hat,tau/2);
        
        % lie splitting

        w_hat=nls_quadratic_nonlinear_part(w_hat,epsilon,tn,tau);
        w_hat=expilaplacian(w_hat,tau);
        tn=tn+tau;
    end

    error_sym_res(j)=norm((u_hat-u_hat_ref).*abs(indexvec).^s);
    error_res(j)=norm((v_hat-u_hat_ref).*abs(indexvec).^s);
    error_lie(j)=norm((w_hat-u_hat_ref).*abs(indexvec).^s);
    error_strang(j)=norm((z_hat-u_hat_ref).*abs(indexvec).^s);
    tau_j(j)=tau;
    end
    
    save(strcat('data/long_time_function_of_tau_epsilon_0-1_H1_data_N_',num2str(N),'.mat'))
end

save(strcat('data/long_time_function_of_tau_epsilon_0-1_H1_data_N_',num2str(N),'.mat'))