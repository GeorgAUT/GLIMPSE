%% Compare the convergence properties of various methods for the SM equation
% Reference values can be computed using compute_reference_solution.m

clear all

%% Choose between two initial conditions of different regularity
%addpath('initial_conditions_smooth/')
addpath('initial_conditions_rough/')

%% Parameters

N=256; % Number of Fourier modes
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
%% Load reference solution
% GMRES Parameters
tol = 1e-13;
maxit = 20;

load('data/reference_solution_rough_N-256_T_1-0.mat','T_ref','h','Nref');
aux=fftpi_mat(T_ref);
T_ref=ifftpi_mat(aux(:,:,Nref/2-N/2+1:Nref/2+N/2))
%% Now compare various time stepping methods
Lmax=5;

hll=zeros(Lmax,1);
err_be=zeros(Lmax,1);
err_chag=zeros(Lmax,1);
err_interpolatory_magnus_O2=zeros(Lmax,1);
err_interpolatory_magnus_O4=zeros(Lmax,1);
err_low_reg=zeros(Lmax,1);

time_be=zeros(Lmax,1);
time_chag=zeros(Lmax,1);
time_interpolatory_magnus_O2=zeros(Lmax,1);
time_interpolatory_magnus_O4=zeros(Lmax,1);
time_low_reg=zeros(Lmax,1);

for l=1:Lmax
    l
    M=l^2;
    h=T/M;
    tol = max(1e-14,h^3);
    maxit = 100;
        
    %  Xie et al. semi-implicit order 1
    Tn=T0;
    tic
    for m=1:M
        % The numerical scheme
        Tn=SME_backward_differ_semi_implicit_first_order(Tn,h,tol,maxit);
    end
    time_be(l)=toc;
    
    T_be=Tn;

    % Xie et al. semi-implicit order 2
    % Initialise scheme
    Tn=T0;
    Tnp1=T0;
    for ll=1:floor(sqrt(M))
        Tnp1=SME_backward_differ_semi_implicit_first_order(Tnp1,h/floor(sqrt(M)),tol,maxit);
    end
    tic
    for m=1:M
        Tnp1_old=Tnp1;
        Tnp1=SME_backward_differ_semi_implicit_second_order(Tn,Tnp1,h,tol,maxit);
        Tn=Tnp1_old;
    end
    time_chag(l)=toc;
    T_chang=Tn;
    
    % Second order interpolatory Magnus integrator
    
    % Get all initial data needed for our scheme
    u0_hat=fftpi(u0);
    u_hat=u0_hat;
    y0=zeros(3,3,N);
    y0(1,:,:)=T0;
    y0(2,:,:)=e10;
    y0(3,:,:)=e20;
    y=y0;
    tic
    for m=1:M

        % Take a timestep with strang
        % Strang
        u_hat_old=u_hat;
        u_hat=strang_NLS(u_hat,h); 
        
        y=interpolatory_magnus_integrator_O2(y,u_hat_old,u_hat,h);
        
    end
    time_interpolatory_magnus_O2(l)=toc;
    T_interpolatory_magnus_O2=y(1,:,:);
    
    % Third order interpolatory Magnus
    % Get all initial data needed for our scheme
    u0_hat=fftpi(u0);
    u_hat=u0_hat;
    y0=zeros(3,3,N);
    y0(1,:,:)=T0;
    y0(2,:,:)=e10;
    y0(3,:,:)=e20;
    y=y0;
    u_hat=u0_hat;
    t1=(1-1/sqrt(3))/2;
    t2=(1+1/sqrt(3))/2;
    tic
    for m=1:M

        % Take a timestep with fourth order NLS splitting method
        u_hat_t1=fourth_order_splitting_NLS(u_hat,h*t1);
        u_hat_t2=fourth_order_splitting_NLS(u_hat_t1,h*(t2-t1));
        u_hat=fourth_order_splitting_NLS(u_hat_t2,h*(1-t2)); 
        
        y=interpolatory_magnus_integrator_O4(y,u_hat_t1,u_hat_t2,h);
        
    end
    time_interpolatory_magnus_O4(l)=toc;
    T_interpolatory_magnus_O4=y(1,:,:);

    % Low regularity integrator

    % Get all initial data needed for our scheme
    u0_hat=fftpi(u0);
    u_hat=u0_hat;
    y0=zeros(3,3,N);
    y0(1,:,:)=T0;
    y0(2,:,:)=e10;
    y0(3,:,:)=e20;
    y=y0;
    tic
    for m=1:M
        M-m
        % Take a timestep with low-regularity integrator
        u_hat_old=u_hat;
        u_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(u_hat,-1/2,h/2),-1/2,h/2);
        
        y=SME_low_reg_integrator(y,u_hat_old,u_hat,h);
        
    end
    time_low_reg(l)=toc;
    T_low_reg=y(1,:,:);
    
    % Now measure L2 error

    hll(l)=h;
    err_be(l)=sqrt(sum(sum(abs(fftpi_mat(T_ref-T_be)).^2)));
    err_chag(l)=sqrt(sum(sum(abs(fftpi_mat(T_ref-T_chang)).^2)));
    err_interpolatory_magnus_O2(l)=sqrt(sum(sum(abs(fftpi_mat(T_ref-T_interpolatory_magnus_O2)).^2)));
    err_interpolatory_magnus_O4(l)=sqrt(sum(sum(abs(fftpi_mat(T_ref-T_interpolatory_magnus_O4)).^2)));
    err_low_reg(l)=sqrt(sum(sum(abs(fftpi_mat(T_ref-T_low_reg)).^2)));

    save(strcat('data/convergence_results_rough_T-',num2str(T),'_N_',num2str(N),'_v4.mat'))
end