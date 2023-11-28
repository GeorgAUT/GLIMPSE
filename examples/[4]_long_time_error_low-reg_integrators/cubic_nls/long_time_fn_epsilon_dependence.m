%% In this example we evaluate the performance of our methods as a function of epsilon
% The results can be plotted using plotting_fn_epsilon.m
%
clear all
rng(10) % Fix random seed for reproducability
%% Parameters

N_test=2^10; % Number of Fourier modes
N_ref=2^12; % Number of Fourier modes in reference solution

N=N_ref;
Xvec=pi*(-N/2+1:N/2)'/N*2;
indexvec=(-N/2+1:N/2)';
indexvec(N/2)=1;

%% Initial conditions of O(1)
% Low-regularity initial conditions
theta=2.0; %regularity parameter for initial condition
u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
u0_hat=u0_hat./indexvec.^theta;
%u0_hat(N/2)=0;

% Smooth initial conditions
% u0=@(x) cos(x)./(2+sin(x));
% u0_eval=u0(Xvec);
% u0_hat=fftpi(u0_eval);

% Equation we consider is i\partial_t u= -\Delta u +\epsilon^2 |u|^2u
%% tau=0.5
Jmax=5;
error_sym_res=zeros(1,Jmax);
error_res=zeros(1,Jmax);
error_strang=zeros(1,Jmax);
error_lie=zeros(1,Jmax);
error_res_long_time=zeros(1,Jmax);
error_res_long_time_second_order=zeros(1,Jmax);
error_res_long_time_second_order_yue=zeros(1,Jmax);
epsilon_j=zeros(1,Jmax); % For collecting timestep sizes
switches=zeros(1,Jmax);
%% Loading saved results

a=input('Would you like to continue the computation that has already started?')

if a==1
    load(strcat('data/long_time_function_of_epsilon_tau_0_05_H2_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'_v1.mat'))
end

s=1 % Sobolev norm of error
for j=1:Jmax

    if switches(j)==0
        %% Compute time-stepping for decreasing epsilon
        Jmax-j
        epsilon=2^(-j+1);
        T=1/epsilon^2; %Final time
        
        %% Compute reference solution with small time-step & resonance-based method
        tau=1e-3;
        M=floor(T/tau);
        tau=T/M;
        u_hat=u0_hat; % Start time-integration from initial condition
        for m=1:M
            M-m
            u_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(u_hat,epsilon^2,tau/2),epsilon^2,tau/2);
        end
        
        u_hat_ref=u_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Save reference solution
        
        
        tau=0.05
        M=floor(T/tau);
        tau=T/M;
        u_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with symmetric resonance-based
        v_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with standard resonance-based
        w_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with strang splitting
        zz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with lie splitting
        zzz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with new resonance_longtime
        zzzz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with new resonance_symmetric second order
        zzzzz_hat=u0_hat(N_ref/2-N_test/2+1:N_ref/2+N_test/2); % Sln computed with new resonance_symmetric second order yue
        tn=0;
        for m=1:M
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
        
        
        N=N_test;
        Xvec=pi*(-N/2+1:N/2)'/N*2;
        indexvec=(-N/2+1:N/2)';
        indexvec(N/2)=1;
        
        
        error_sym_res(j)=norm((u_hat-u_hat_ref).*abs(indexvec).^s);
        error_res(j)=norm((v_hat-u_hat_ref).*abs(indexvec).^s);
        error_strang(j)=norm((w_hat-u_hat_ref).*abs(indexvec).^s);
        error_lie(j)=norm((zz_hat-u_hat_ref).*abs(indexvec).^s);
        error_res_long_time(j)=norm((zzz_hat-u_hat_ref).*abs(indexvec).^s);
        error_res_long_time_second_order(j)=norm((zzzz_hat-u_hat_ref).*abs(indexvec).^s);
        error_res_long_time_second_order_yue(j)=norm((zzzzz_hat-u_hat_ref).*abs(indexvec).^s);
        
        epsilon_j(j)=epsilon;
        
        switches(j)=1;
    end

    save(strcat('data/long_time_function_of_epsilon_tau_0_05_H2_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'_v1.mat'))
end

save(strcat('data/long_time_function_of_epsilon_tau_0_05_H2_data_Ntest_',num2str(N_test),'_Nref_',num2str(N_ref),'_v1.mat'))

