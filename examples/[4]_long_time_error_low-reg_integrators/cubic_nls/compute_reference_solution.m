%% Precomputation of reference solution for convergence tests

clear all
rng(10) % Fix random seed for reproducability
%% Parameters

N=2^8; % Number of Fourier modes
Xvec=pi*(-N/2+1:N/2)'/N*2;
indexvec=(-N/2+1:N/2)';
indexvec(N/2)=1;
%% Initial conditions of O(1)

% Low-regularity initial conditions
theta=2.5; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;

% Smooth initial conditions
% u0=@(x) cos(x)./(2+sin(x));
% u0_eval=u0(Xvec);
% u0_hat=fftpi(u0_eval);

%% Equation we consider is i\partial_t u= -\Delta u +\epsilon^2 |u|^2u
% Epsilon 0.1

epsilon=0.5;
T=1/epsilon^2; %Final time

%% Compute reference solution with small time-step & resonance-based method
tau=0.2*1e-3;
M=floor(T/tau);
tau=T/M;
u_hat=u0_hat; % Start time-integration from initial condition
for m=1:M
    M-m
    u_hat=cubic_nls_resonance_based_first_order_NRSLI1_adjoint(cubic_nls_resonance_based_first_order_NRSLI1(u_hat,epsilon^2,tau/2),epsilon^2,tau/2);
    %input('')
end

u_hat_ref=u_hat; % Save reference solution
N_ref=N;

save(strcat('data/reference_solution_epsilon_0-5_H2_N_256.mat'),'u0_hat','u_hat_ref','N_ref')
