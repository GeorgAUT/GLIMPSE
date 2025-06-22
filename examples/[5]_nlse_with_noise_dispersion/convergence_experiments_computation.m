%% Here we evaluate  the pathwise error of our integrator against prior work
% NLSE equation with white noise dispersion
%
% The code allows for multiple sets of parameters to be specified and run
% in sequence. The data is saved and can then be plotted using
% plot_convergence_results.m
%
clear all
rng(10) % Fix random seed for reproducability
%% add the correct paths in Matlab to all the submodules

run('../../src/set_paths.m')

%% Operational parameter vector
exp_settings={{'H4',2^6,2^7,1.0,2^6,2000,0.1,1}};% Format: {'text_reg',N,Nref,mu,N_obs,Ntau,epsilon,X}

for jjj=1:size(exp_settings,2)
    text_reg=exp_settings{jjj}{1}; % Regularity parameter: options = 'Cinf', 'H2', 'H3', 'H4' 
    N=exp_settings{jjj}{2}; % Number of Fourier modes
    Nref=exp_settings{jjj}{3}; % Number of Fourier modes in reference solution
    mu=exp_settings{jjj}{4}; % Equation we consider is i\partial_t u= -\Delta u +mu|u|^2u
    N_obs=exp_settings{jjj}{5}; % Number of observations in Wong-Zakai, i.e. \delta=T/N_obs
    Ntau=exp_settings{jjj}{6}; % Number of different timesteps used to plot convergence graphs

    epsilon=exp_settings{jjj}{7}; % Scaling of L^2 norm of initial data
    X=exp_settings{jjj}{8}; % Number of simulation paths
    
    
    Xvec=pi*(-N/2+1:N/2)'/N*2; % Spatial mesh
    indexvec=(-N/2+1:N/2)'; % Vector of Fourier indices
    indexvec(N/2)=1; % Correct zero index so can use to rescale Fourier coefficients
    
    Xvecref=pi*(-Nref/2+1:Nref/2)'/Nref*2; % Spatial mesh
    indexvecref=(-Nref/2+1:Nref/2)'; % Vector of Fourier indices
    indexvecref(Nref/2)=1; % Correct zero index so can use to rescale Fourier coefficients
    
    T=1.0; % Final time
    
    s=1; % Sobolev norm for scaling of initial data and for measuring convergence
    
    %% Initial data
    switch text_reg
        case 'Cinf'
            u0_hat_ref=fftpi(sin(Xvecref)./(2+cos(Xvecref))); %smooth initial data
        case 'H1'
            theta=1.0; %regularity parameter for initial condition
            u0_hat_ref=rand(size(Xvecref))+i*rand(size(Xvecref));
            u0_hat_ref=u0_hat_ref./indexvecref.^theta;
        case 'H2'
            theta=2.0; %regularity parameter for initial condition
            u0_hat_ref=rand(size(Xvecref))+i*rand(size(Xvecref));
            u0_hat_ref=u0_hat_ref./indexvecref.^theta;
        case 'H3'
            theta=3.0; %regularity parameter for initial condition
            u0_hat_ref=rand(size(Xvecref))+i*rand(size(Xvecref));
            u0_hat_ref=u0_hat_ref./indexvecref.^theta;
        case 'H4'
            theta=4.0; %regularity parameter for initial condition
            u0_hat_ref=rand(size(Xvecref))+i*rand(size(Xvecref));
            u0_hat_ref=u0_hat_ref./indexvecref.^theta;
    end
    u0_hat_ref(Nref/2)=rand(1);
    u0_hat_ref=u0_hat_ref/norm(u0_hat_ref.*indexvecref.^s)*epsilon; % normalisation of initial conditions
    
    %% Compute error for multiple time steps
        
    H1_low_reg_new=zeros(1,Ntau);
    H1_low_reg_symp=zeros(1,Ntau);
    H1_low_reg_hybrid=zeros(1,Ntau);
    H1_splitting=zeros(1,Ntau);
    H1_relaxedCN=zeros(1,Ntau);
    H1_exponential=zeros(1,Ntau);
    H1_MC=zeros(1,Ntau);
    
    tauj=zeros(Ntau,1);
    mmm_current=0;
    for mmm=1:X
        mmm % Counter for number of samples
    
        % Observe Brownian motion
        deltat_obs=T/N_obs;
        w=sqrt(deltat_obs)*randn(1,N_obs);
        B=zeros(1,N_obs+1);
        
        for j=1:N_obs
            j;
            B(:,j+1)=sum(w(1:j),2);
        end
        
        % Define linear interpolant
        Bfunc=@(x) interp1(linspace(0,T,N_obs+1),B,x);
    
        x_plot=linspace(0,T,100);
        
        figure(1)
        plot(x_plot,Bfunc(x_plot))
        
        %% Compute reference solution
        Mref=2^15
        tau_ref=T/Mref;
        
        u_hat_ref=u0_hat_ref;
        
        recompute=1; % Always recompute the reference solution
        filename=strcat('data/reference_solution_',text_reg,'_Nref_',num2str(Nref),'_T_',num2str(T),'_Nobs_',num2str(N_obs),'_Mref_',num2str(Mref),'.mat');
        
        if isfile(filename) & 1-recompute
            load(filename,'u_hat_ref','u0_hat','Bfunc')
        else
            u_hat_ref=u0_hat_ref;
            for m=1:Mref
                Mref-m
                % Lie splitting:
                psin=[0,Bfunc(m*tau_ref)-Bfunc((m-1)*tau_ref)];
                u_hat_ref=expilaplacian(u_hat_ref,psin(2));
                u_hat_ref=nls_nonlinear_part(u_hat_ref,-mu*tau_ref);
            end
            %save(filename)
        end
    
        %% Downsampling
        u0_hat=u0_hat_ref(Nref/2-N/2+1:Nref/2+N/2);
        u_hat_ref=u_hat_ref(Nref/2-N/2+1:Nref/2+N/2);
        indexvec=indexvecref(Nref/2-N/2+1:Nref/2+N/2);
        
        %% Now try for a few different parameters
        
        parfor jj=1:Ntau
            Ntau-jj % Output where in the calculation we are
    
            M=jj; % Define timestep (through M = total number of time steps)
            M=T*M;
            tau=T/M;
            tauj(jj)=tau;
        
            u_hat=u0_hat; % Initialise Lie splitting
            v_hat=u0_hat; % Initialise Exponential integrator
            vv_hat=u0_hat; % Initialise Exponential integrator
            znew_hat=u0_hat; % Low-regularity integrator SDLR-I
            zzz_hat=u0_hat; % Initialise CN
    
            % Auxiliary computation:
    
            N1=N_obs; % Fix Wong-Zakai a priori
            tau=T/M;
            tauvec=tau*(0:M);
            varepsilon=T/N_obs;
            varepsilonvec=varepsilon*(0:N_obs);
            
            timeintervalscombined=[tauvec,varepsilonvec];
    
            timeintervalscombined=sort(timeintervalscombined);
            repeats = diff(timeintervalscombined) <1e-14; 
            timeintervalscombined(repeats) = [];
    
            zzz_hat_mtauhalf=expilaplacian(zzz_hat,Bfunc(tau/2));
            zzz_hat_mtauhalf=nls_nonlinear_part(zzz_hat_mtauhalf,-mu*tau/2);
            Gamma_hat=fftpi(abs(ifftpi(zzz_hat_mtauhalf)).^2);% Initialise aux variable in relaxed CN
    
            % Compute time steps
    
            for m=1:M
                psin=[0,Bfunc(m*tau)-Bfunc((m-1)*tau)];
    
                % Splitting
                u_hat=expilaplacian(u_hat,psin(2));
                u_hat=nls_nonlinear_part(u_hat,-mu*tau);
        
                % relaxed CN
                [zzz_hat,Gamma_hat]=cubic_nls_crank_nicolson_relaxed_variable_dispersion(zzz_hat,Gamma_hat,-mu,tau,psin(2));
    
                
    
            end
    
            for m=1:M
    
                tn=tau*(m-1);
                tnp1=tau*m;
    
                % Find timeintervals inside [tn,tnp1]
                index_list=find((timeintervalscombined>=tn-1e-15).*(timeintervalscombined<=tnp1+1e-15));
    
                local_t_list=timeintervalscombined(index_list); % List of all entries in that time interval
                psin=Bfunc(local_t_list)-Bfunc(tn);
    
                % First low-reg integrator
                znew_hat=cubic_nls_resonance_based_first_order_variabledispersion(znew_hat,-mu,psin,local_t_list);
                
                % Exponential integrator new
                taulist=local_t_list(2:end)-local_t_list(1:end-1);
                psin=Bfunc(local_t_list)-Bfunc(tn);
                vv_phys=ifftpi(vv_hat);
                aux1=fftpi(mu*abs(vv_phys).^2.*vv_phys);
                vv_hat=expilaplacian(vv_hat,psin(end))+...
                    i*expilaplacian(sum(taulist.*expilaplacian(phi1laplacian(aux1,-(psin(2:end)-psin(1:end-1))),-psin(1:end-1)),2),psin(end));
            
            end
    
            % Evaluate the pathwise error
           
            H1_low_reg_new(:,jj)=H1_low_reg_new(:,jj)+(sum((abs(znew_hat-u_hat_ref).*indexvec.^s).^2,1));
            H1_splitting(:,jj)=H1_splitting(:,jj)+(sum((abs(u_hat-u_hat_ref).*indexvec.^s).^2,1));
            H1_relaxedCN(:,jj)=H1_relaxedCN(:,jj)+(sum((abs(zzz_hat-u_hat_ref).*indexvec.^s).^2,1));
            H1_exponential(:,jj)=H1_exponential(:,jj)+(sum((abs(vv_hat-u_hat_ref).*indexvec.^s).^2,1));
        end
        %% Save data
        save(strcat('data/convergence_plot_NLSE_Tend_',num2str(T),'_',text_reg,'_N_',num2str(N),'_Nref_',num2str(Nref),'_Ntau_',num2str(Ntau),'_Nobs_',num2str(N_obs),'_epsilon_',strrep(num2str(epsilon),'.','-'),'_X_',num2str(X),'_piecewise_linear_varorig.mat'))
        % Saving here allows for more dynamic averaging
    end
end