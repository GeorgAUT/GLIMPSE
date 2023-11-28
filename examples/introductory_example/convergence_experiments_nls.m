%% This is a simple example to demonstrate how the methods can be used

clear all
rng(1) % Fix random seed for reproducability

%% Parameters
%
N=2048; % Number of Fourier modes
Xvec=pi*(-N/2+1:N/2)'/N*2; % Spacial mesh
indexvec=(-N/2+1:N/2)'; % Vector of Fourier indices
indexvec(N/2)=1; % Correct zero index so can use to rescale Fourier coefficients
mu=1.0; % Equation we consider is i\partial_t u= -\Delta u +|u|^2u
T=1.0; % Final time

%% Initial conditions (H^1 data)
%
theta=1.0; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;

u0_hat=u0_hat/norm(u0_hat); % normalisation
%% Reference solution computed with second order scheme from Bruned & Schratz
Mref=1e5;
tauref=T/Mref;

recompute=input('\n If you would like to recompute the reference sln please enter 1: ');

if recompute==1
    u_hat_ref=u0_hat;
    for m=1:Mref
        Mref-m
        
        % Standard second order
        u_hat_ref=cubic_nls_resonance_based_second_order_BS22(u_hat_ref,tauref);
    
    end
    save('data/reference_solution.mat','u_hat_ref');
else
    load('data/reference_solution.mat','u_hat_ref');
end

%% Compute L^2 error in time-stepping for decreasing tau
Jmax=100;

error_first_order=zeros(1,Jmax);
error_brunedschratz=zeros(1,Jmax);
error_symmetric=zeros(1,Jmax);
error_strang=zeros(1,Jmax);
error_lie=zeros(1,Jmax);

cputime_first_order=zeros(1,Jmax);
cputime_brunedschratz=zeros(1,Jmax);
cputime_symmetric=zeros(1,Jmax);
cputime_strang=zeros(1,Jmax);
cputime_lie=zeros(1,Jmax);
tau_j=zeros(1,Jmax);

parfor j=1:Jmax
    Jmax-j
    tau=1/(10*j)^1.5;
    M=floor(T/tau);
    tau=T/M;

    u_hat=u0_hat;
    v_hat=u0_hat;
    w_hat=u0_hat;
    z_hat=u0_hat;
    zz_hat=u0_hat;
    
    % Standard first order
    tic
    for m=1:M
        u_hat=cubic_nls_resonance_based_first_order(u_hat,mu,tau);
    end
    cputime_first_order(j)=toc; % CPU-time
    error_first_order(j)=norm(u_hat-u_hat_ref); % L^2-error

    % Symmetrised second order
    tic
    for m=1:M
        v_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(v_hat,mu,tau/2),mu,tau/2);
    end
    cputime_symmetric(j)=toc; % CPU-time
    error_symmetric(j)=norm(v_hat-u_hat_ref); % L^2-error
        
    % Standard second order
    tic
    for m=1:M
        w_hat=cubic_nls_resonance_based_second_order_BS22(w_hat,tau);
    end
    cputime_brunedschratz(j)=toc; % CPU-time
    error_brunedschratz(j)=norm(w_hat-u_hat_ref); % L^2-error
    
    % Lie
    tic
    for m=1:M
        z_hat=expilaplacian(nls_nonlinear_part(z_hat,tau),tau);
    end
    cputime_lie(j)=toc; % CPU-time
    error_lie(j)=norm(z_hat-u_hat_ref); % L^2-error

    % Strang
    tic
    for m=1:M
        zz_hat=expilaplacian(nls_nonlinear_part(expilaplacian(zz_hat,tau/2),tau),tau/2);
    end
    cputime_strang(j)=toc;
    error_strang(j)=norm(zz_hat-u_hat_ref);

    tau_j(j)=tau; % Store time step corresponding to experiment

end    


save('data/numerical_experiments_H1_data.mat')