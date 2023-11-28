%% Compute reference solution to SM equation for convergence experiments

clear all

%% Choose between two initial conditions of different regularity
%addpath('initial_conditions_smooth/')
addpath('initial_conditions_rough/')

%% Parameters

N=2^8; % Number of Fourier modes
Xvec=pi*(-N/2+1:N/2)/N*2;
T=1.0;

%% Set up the evolution vectors/matrix...

T0=zeros(1,3,N);
e10=zeros(1,3,N);
e20=zeros(1,3,N);

for j=1:N
    T0(:,:,j)=T0_cts(Xvec(j));
    e10(:,:,j)=e10_cts(Xvec(j));
    e20(:,:,j)=e20_cts(Xvec(j));
end

u0=zeros(size(Xvec'));
for j=1:N
    u0(j)=u0_cts(Xvec(j));
end
%% Compute reference solution
M=2^8; % Number of time steps
h=T/M;
% GMRES Parameters
tol = 1e-13;
maxit = 20;

%% Low-regularity_integrator
% Get all initial data needed for our scheme
u0_hat=fftpi(u0);
u_hat=u0_hat;
y0=zeros(3,3,N);
y0(1,:,:)=T0;
y0(2,:,:)=e10;
y0(3,:,:)=e20;
y=y0;
u_hat=u0_hat;
   for m=1:M
        M-m
        % Take a timestep with symmetric low-regularity integrator
        u_hat_old=u_hat;
        u_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(u_hat,-1/2,h/2),-1/2,h/2);
        
        % Take a step with FLowRH transform
        y=SME_low_reg_integrator(y,u_hat_old,u_hat,h);
        
    end
T_ref=y(1,:,:);
Nref=N;

save(strcat('data/reference_solution_rough_N-',num2str(N),'_T_1-0.mat'),'T_ref','h','Nref');