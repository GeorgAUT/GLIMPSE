%% Comparison of error over long time periods for fixed tau
% To plot the results please use plotting_time_evolution_error.m
%
clear all
addpath('modules/')
rng(10) % Fix random see for reproducability
%% Parameters

N_test=2^13; % Number of Fourier modes
N_ref=2^14;

N=N_ref;

Xvec=pi*(-N/2+1:N/2)'/N*2;
indexvec=(-N/2+1:N/2)';
indexvec(N/2)=1;


%% Initial conditions of O(1)
% Low-regularity initial conditions
theta=2.0; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;
% u0_hat(N/2)=0;
% Smooth initial conditions
% u0=@(x) cos(x)./(2+sin(x));
% u0_eval=u0(Xvec);
% u0_hat=fftpi(u0_eval);

%% Equation we consider is i\partial_t u= -\Delta u +\epsilon^2 |u|^2u
% Epsilon=0.1
epsilon=0.1;
T=1/epsilon^2; %Final time

%% Computing the evolution of the error

% Parameters
tau=1e-2;
M=floor(T/tau);
tau=T/M;

% Arrays storing error values and time values
t_vec=zeros(1,M);
error_sym_res=zeros(1,M);
error_res=zeros(1,M);
error_lie=zeros(1,M);
error_strang=zeros(1,M);
error_res_long_time=zeros(1,M);
error_res_long_time_second_order=zeros(1,M);

% Additional refinement of temporal grid for reference solution
ref_extra=100;
tau_ref=tau/ref_extra;

% Number of spatial Fourier modes in the tests
N=N_test;
indexvec=(-N/2+1:N/2)';
indexvec(N/2)=1;

%% Initiation and time propagation
u_ref_hat=u0_hat; % reference solution values

u_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with symmetric resonance-based
v_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with standard resonance-based
w_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with strang splitting
zz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with lie splitting
zzz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with new resonance_longtime
zzzz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with new resonance_longtime second order
zzzzz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with new resonance_longtime second order

s=1.0;

for m=1:M
    M-m
    for mm=1:ref_extra
        u_ref_hat=cubic_nls_resonance_based_first_order_NRSLI1_adjoint(cubic_nls_resonance_based_first_order_NRSLI1(u_ref_hat,epsilon^2,tau/2/ref_extra),epsilon^2,tau/2/ref_extra);
    end
    tn=0;
    % Resonance-based methods 
    u_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(u_hat,epsilon^2,tau/2),epsilon^2,tau/2);
    v_hat=cubic_nls_resonance_based_first_order(v_hat,epsilon^2,tau);
    zzz_hat=cubic_nls_resonance_based_first_order_NRSLI1(zzz_hat,epsilon^2,tau);
    zzzzz_hat=cubic_nls_resonance_based_second_order_NRSLI2(zzzzz_hat,epsilon^2,tau);
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
    tn=tn+tau;

    % Evaluate the errors
    error_sym_res(m)=norm((u_hat-u_ref_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2)).*abs(indexvec).^s);
    error_res(m)=norm((v_hat-u_ref_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2)).*abs(indexvec).^s);
    error_strang(m)=norm((w_hat-u_ref_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2)).*abs(indexvec).^s);
    error_lie(m)=norm((zz_hat-u_ref_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2)).*abs(indexvec).^s);
    error_res_long_time(m)=norm((zzz_hat-u_ref_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2)).*abs(indexvec).^s);
    error_res_long_time_second_order(m)=norm((zzzz_hat-u_ref_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2)).*abs(indexvec).^s);
    t_vec(m)=m*tau;

    if mod(m,10)==0 % store the results
        save(strcat('data/error_evolution_eps_0-1_H2_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'.mat'))
    end
end

save(strcat('data/error_evolution_eps_0-1_H2_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'.mat'))