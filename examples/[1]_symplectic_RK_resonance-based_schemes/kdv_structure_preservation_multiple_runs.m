%% Here we test the structure preservation of our symplectic method against prior work
% KdV equation
% The code allows for multiple sets of parameters to be specified and run
% in sequence. The data is saved and can then be plotted using
% kdv_plot_structure_preservation.m
clear all
rng(25) % Fix random seed for reproducability 10

%% Operational vector
exp_settings={{'Cinf',64,400,0.1},...
    {'Cinf',64,40,0.02}}; % Format: {'text_reg',N,T,tau}


nonlinear_part=@(u) 1/2*dx(conv1(u,u)); % Nonlinear part for Strang splitting


for jjj=1:size(exp_settings,2)

    text_reg=exp_settings{jjj}{1};
    N=exp_settings{jjj}{2};
    T=exp_settings{jjj}{3};
    tau=exp_settings{jjj}{4};

    
    %% Parameters
    %
    Xvec=pi*(-N/2+1:N/2)'/N*2; % Spacial mesh
    indexvec=(-N/2+1:N/2)'; % Vector of Fourier indices
    indexvec(N/2)=1; % Correct zero index so can use to rescale Fourier coefficients
    mu=1.0; % Equation we consider is i\partial_t u= -\Delta u +|u|^2u
    
    
    %% Initial conditions (H^1 data)
    %
    switch text_reg
        case 'Cinf'
            u0_hat=fftpi(sin(Xvec)./(2+cos(Xvec))); %smooth initial data
        case 'H2'
            theta=2.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./indexvec.^theta;
        case 'H3'
            theta=3.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./indexvec.^theta;
        case 'H4'
            theta=4.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./indexvec.^theta;
    end
    
    u0_hat=u0_hat/norm(u0_hat)/10; % normalisation
    
    % Define Hamiltonian
    Ham=@(zzz_hat) sum(ifftpi(zzz_hat).^3+3*ifftpi(dx(zzz_hat)).^2)/N;
    %% Compute the normalisation & Hamiltonian after each step
    M=floor(T/tau);
    tau=T/M;
    
    momentum_brunedschratz=zeros(1,M);
    momentum_strang=zeros(1,M);
    momentum_lawson=zeros(1,M);
    momentum_res_midpoint=zeros(1,M);
    hamiltonian_brunedschratz=zeros(1,M);
    hamiltonian_strang=zeros(1,M);
    hamiltonian_lawson=zeros(1,M);
    hamiltonian_res_midpoint=zeros(1,M);
    t_j=zeros(1,M);
    
    u_hat=u0_hat;
    v_hat=u0_hat;
    w_hat=u0_hat;
    ww_hat=u0_hat;
    www_hat=u0_hat;
    wwww_hat=u0_hat;
    z_hat=u0_hat;
    zz_hat=u0_hat;
    zzz_hat=u0_hat;
    
    for j=1:M
        M-j
    
        % Res-based midpoint
        wwww_hat=kdv_resonance_based_midpoint(wwww_hat,tau);
    
        momentum_res_midpoint(j)=(sum(ifftpi(wwww_hat).^2))/N;
        hamiltonian_res_midpoint(j)=Ham(wwww_hat);
    
        % Standard second order
        w_hat=kdv_resonance_based_second_order_bruned_schratz(w_hat,tau);
    
        momentum_brunedschratz(j)=(sum(ifftpi(w_hat).^2))/N;
        hamiltonian_brunedschratz(j)=Ham(w_hat);
        
        % Strang
        zz_hat=expmdx3(zz_hat,tau/2);
        hnl=tau/1000;
        
        for mm=1:1000
            k1=nonlinear_part(zz_hat);
            k2=nonlinear_part(zz_hat+hnl*k1/2);
            k3=nonlinear_part(zz_hat+hnl*k2/2);
            k4=nonlinear_part(zz_hat+hnl*k3);
            zz_hat=zz_hat+1/6*hnl*(k1+2*k2+2*k3+k4);
        end
        
        zz_hat=expmdx3(zz_hat,tau/2);
    
        momentum_strang(j)=(sum(ifftpi(zz_hat).^2))/N;
        hamiltonian_strang(j)=Ham(zz_hat);
    
        % Lawson
        zzz_hat=kdv_onestageLawson1(zzz_hat,tau);
    
        momentum_lawson(j)=(sum(ifftpi(zzz_hat).^2))/N;
        hamiltonian_lawson(j)=Ham(zzz_hat);
    
        t_j(j)=j*tau; % Store time step corresponding to experiment

        if mod(j,1000)==0
            save(strcat('data/kdv_structure_preservation_scaled_',text_reg,'_N_',num2str(N),'_T_',num2str(T),'_tau_',strrep(num2str(tau),'.','-'),'.mat'))
        end
    end
    
    
    save(strcat('data/kdv_structure_preservation_scaled_',text_reg,'_N_',num2str(N),'_T_',num2str(T),'_tau_',strrep(num2str(tau),'.','-'),'.mat'))

end

