%% Here we test the convergence of our symmetric methods against prior work
% NLS equation
% The code allows for multiple sets of parameters to be specified and run
% in sequence. The data is saved and can then be plotted using
% plotting_convergence_graphs.m

clear all
rng(10) % Fix random seed for reproducability
%% Operational vector
exp_settings={{'H2',2^8,2^6,10},...
    {'Cinf',2^8,2^6,10}}; % Format: {'text_reg',Nref,N,Jmax}

for jjj=1:size(exp_settings,2)

    text_reg=exp_settings{jjj}{1};
    Nref=exp_settings{jjj}{2};
    N=exp_settings{jjj}{3};
    Jmax=exp_settings{jjj}{4};
    %% Parameters
    %
    Xvec=pi*(-Nref/2+1:Nref/2)'/Nref*2; % Spacial mesh
    indexvec=(-Nref/2+1:Nref/2)'; % Vector of Fourier indices
    indexvec(Nref/2)=1; % Correct zero index so can use to rescale Fourier coefficients
    mu=1.0; % Equation we consider is i\partial_t u= -\Delta u +|u|^2u
    T=1.0; % Final time


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
    
    u0_hat=u0_hat/norm(u0_hat); % normalisation
    
    %% Reference solution computed with second order scheme from Bruned & Schratz
    Mref=1e4;
    tauref=T/Mref;
    
    recompute=1;%input('\n If you would like to recompute the reference sln please enter 1: ');
    filename=strcat('data/reference_solution_',text_reg,'_N_',num2str(Nref),'_T_',num2str(T),'.mat');
    
    if isfile(filename)
        load(filename,'u_hat_ref','u0_hat')
    else
        u_hat_ref=u0_hat;
        for m=1:Mref
            Mref-m
            
            % Standard second order
            u_hat_ref=cubic_nls_resonance_based_second_order_BS22(u_hat_ref,tauref);
        
        end
        save(filename)
    end
    
    %% Downsampling
    u0_hat=u0_hat(Nref/2-N/2+1:Nref/2+N/2);
    u_hat_ref=u_hat_ref(Nref/2-N/2+1:Nref/2+N/2);
    indexvec=indexvec(Nref/2-N/2+1:Nref/2+N/2);

    %% Compute L^2 error in time-stepping for decreasing tau
    
    error_first_order=zeros(1,Jmax);
    error_new_symm_first_order=zeros(1,Jmax);
    error_new_symm_second_order=zeros(1,Jmax);
    error_brunedschratz=zeros(1,Jmax);
    error_symmetric=zeros(1,Jmax);
    error_strang=zeros(1,Jmax);
    error_lie=zeros(1,Jmax);
    
    cputime_first_order=zeros(1,Jmax);
    cputime_new_symm_first_order=zeros(1,Jmax);
    cputime_new_symm_second_order=zeros(1,Jmax);
    cputime_brunedschratz=zeros(1,Jmax);
    cputime_symmetric=zeros(1,Jmax);
    cputime_strang=zeros(1,Jmax);
    cputime_lie=zeros(1,Jmax);
    tau_j=zeros(1,Jmax);

    s=1.0; % Order of Sobolev norm
    
    filename=strcat('data/convergence_data_nlse_',text_reg,'_Nref_',num2str(Nref),'_N_',num2str(N),'_T_',num2str(T),'_Jmax_',num2str(Jmax),'.mat');

    parfor j=1:Jmax
        Jmax-j
        tau=1/(10*j)^1.5;
        M=floor(T/tau);
        tau=T/M;
    
        u_hat=u0_hat;
        v_hat=u0_hat;
        w_hat=u0_hat;
        ww_hat=u0_hat;
        www_hat=u0_hat;
        z_hat=u0_hat;
        zz_hat=u0_hat;
        
        % Standard first order
        tic
        for m=1:M
            u_hat=cubic_nls_resonance_based_first_order(u_hat,mu,tau);
        end
        cputime_first_order(j)=toc; % CPU-time
        error_first_order(j)=norm((u_hat-u_hat_ref).*indexvec.^s); % H^s-error
    
        % New symmetric first order
        tic
        for m=1:M
            ww_hat=cubic_nls_resonance_based_first_order_symm_new(ww_hat,mu,tau);
        end
        cputime_new_symm_first_order(j)=toc; % CPU-time
        error_new_symm_first_order(j)=norm((ww_hat-u_hat_ref).*indexvec.^s); % H^s-error
    
        % New symmetric second order
        tic
        for m=1:M
            www_hat=cubic_nls_resonance_based_second_order_symm_new(www_hat,mu,tau);
        end
        cputime_new_symm_second_order(j)=toc; % CPU-time
        error_new_symm_second_order(j)=norm((www_hat-u_hat_ref).*indexvec.^s); % H^s-error
    
        % Symmetrised second order
        tic
        for m=1:M
            v_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(v_hat,mu,tau/2),mu,tau/2);
        end
        cputime_symmetric(j)=toc; % CPU-time
        error_symmetric(j)=norm((v_hat-u_hat_ref).*indexvec.^s); % H^s-error
            
        % Standard second order
        tic
        for m=1:M
            w_hat=cubic_nls_resonance_based_second_order_BS22(w_hat,tau);
        end
        cputime_brunedschratz(j)=toc; % CPU-time
        error_brunedschratz(j)=norm((w_hat-u_hat_ref).*indexvec.^s); % H^s-error
        
        % Lie
        tic
        for m=1:M
            z_hat=expilaplacian(nls_nonlinear_part(z_hat,tau),tau);
        end
        cputime_lie(j)=toc; % CPU-time
        error_lie(j)=norm((z_hat-u_hat_ref).*indexvec.^s); % H^s-error
    
        % Strang
        tic
        for m=1:M
            zz_hat=expilaplacian(nls_nonlinear_part(expilaplacian(zz_hat,tau/2),tau),tau/2);
        end
        cputime_strang(j)=toc;
        error_strang(j)=norm((zz_hat-u_hat_ref).*indexvec.^s);

        tau_j(j)=tau; % Store time step corresponding to experiment
        
        %save(filename);
    end 

    save(filename);
end