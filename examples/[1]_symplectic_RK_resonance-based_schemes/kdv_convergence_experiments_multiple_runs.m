%% Here we test the convergence of our symplectic method against prior work
% KdV equation
% The code allows for multiple sets of parameters to be specified and run
% in sequence. The data is saved and can then be plotted using
% kdv_plot_convergence_rates.m
%
clear all
rng(10) % Fix random seed for reproducability
%% Operational vector
exp_settings={{'H4',2^10,64,20},...
    {'H3',2^10,128,20}}; % Format: {'text_reg',Nref,N,Jmax}

nonlinear_part=@(u) 1/2*dx(conv1(u,u)); % define function for KdV splitting

for jjj=1:size(exp_settings,2)


    %% Parameters from operational vector
    %
    text_reg=exp_settings{jjj}{1};
    Nref=exp_settings{jjj}{2};
    N=exp_settings{jjj}{3};
    Jmax=exp_settings{jjj}{4};

    Xvec=pi*(-Nref/2+1:Nref/2)'/Nref*2; % Spacial mesh
    indexvec=(-Nref/2+1:Nref/2)'; % Vector of Fourier indices
    indexvec(Nref/2)=1; % Correct zero index so can use to rescale Fourier coefficients
    mu=1.0; % Equation we consider is i\partial_t u= -\Delta u +|u|^2u
    T=1.0; % Final time


    %% Initial conditions
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
        case 'H5'
            theta=5.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./indexvec.^theta;
    end
    u0_hat(Nref/2)=0;
    u0_hat=u0_hat/norm(u0_hat); % normalisation
    
    %% Reference solution computed with second order scheme from Bruned & Schratz
    Mref=1e4;
    tauref=T/Mref;
    
    filename=strcat('data/kdv_reference_solution_',text_reg,'_N_',num2str(Nref),'_T_',num2str(T),'.mat');
    
    if isfile(filename)
        load(filename,'u_hat_ref','u0_hat')
    else
        u_hat_ref=u0_hat;
        for m=1:Mref
            Mref-m
            
            % Standard second order Scheme
            u_hat_ref=kdv_resonance_based_second_order_bruned_schratz(u_hat_ref,tauref);
        
        end
        save(filename)
    end

    %% Downsampling
    u0_hat=u0_hat(Nref/2-N/2+1:Nref/2+N/2);
    u_hat_ref=u_hat_ref(Nref/2-N/2+1:Nref/2+N/2);
    indexvec=indexvec(Nref/2-N/2+1:Nref/2+N/2);


    %% Compute L^2 error in time-stepping for decreasing tau
    
    err_res_midpoint=zeros(Jmax,1);
    err_exp=zeros(Jmax,1);
    err_standard=zeros(Jmax,1);
    err_brunedschratz=zeros(Jmax,1);
    err_lawson=zeros(Jmax,1);
    err_rk4=zeros(Jmax,1);
    err_strang=zeros(Jmax,1);
    
    cputime_res_midpoint=zeros(Jmax,1);
    cputime_exp=zeros(Jmax,1);
    cputime_standard=zeros(Jmax,1);
    cputime_brunedschratz=zeros(Jmax,1);
    cputime_lawson=zeros(Jmax,1);
    hl=zeros(Jmax,1);
    theta=1.0;

    filename=strcat('data/kdv_convergence_data_',text_reg,'_Nref_',num2str(Nref),'_N_',num2str(N),'_T_',num2str(T),'_Jmax_',num2str(Jmax),'.mat');

    parfor l=1:Jmax
        M=10*l;
        Jmax-l
        h=T/M;

        % Schratz and Hofmanova
        
        w_hat=zeros(1,2*M);
        w_hat=u0_hat;
        indexvec=(-N/2+1:N/2);
        indexvec(N/2)=1;
        tic
        for m=1:M
            M-m;
            wstep=w_hat;
            w_hat=kdv_resonance_based_first_order(wstep,h);
    
        end
        cputime_standard(l)=toc;
        err_standard(l)=norm(abs((w_hat-u_hat_ref)).*abs(indexvec).^theta);
        
        % Bruned and Schratz
        
        w_hat=zeros(1,2*M);
        w_hat=u0_hat;
        indexvec=(-N/2+1:N/2);
        indexvec(N/2)=1;
        tic
        for m=1:M
            M-m;
            wstep=w_hat;
            w_hat=kdv_resonance_based_second_order_bruned_schratz(wstep,h);
    
        end
        cputime_brunedschratz(l)=toc;
        err_brunedschratz(l)=norm(abs((w_hat-u_hat_ref)).*abs(indexvec).^theta);
        

        % New resonance-based midpoint
        
        z_hat=zeros(size(w_hat));
        z_hat=u0_hat;
        tic
        for m=1:M
             M-m;
             wstep=z_hat;
             z_hat=kdv_resonance_based_midpoint(wstep,h);
        end
        cputime_res_midpoint(l)=toc;
        err_res_midpoint(l)=norm(abs((z_hat-u_hat_ref)).*abs(indexvec).^theta);
        
        % Strang splitting
        zz_hat=u0_hat;
        for m=1:M
            M-m;
            zz_hat=expmdx3(zz_hat,h/2);
            hnl=h/1000;
            
            for mm=1:1000
                k1=nonlinear_part(zz_hat);
                k2=nonlinear_part(zz_hat+hnl*k1/2);
                k3=nonlinear_part(zz_hat+hnl*k2/2);
                k4=nonlinear_part(zz_hat+hnl*k3);
                zz_hat=zz_hat+1/6*hnl*(k1+2*k2+2*k3+k4);
            end
            
            zz_hat=expmdx3(zz_hat,h/2);
        end
        err_strang(l)=norm(abs((zz_hat-u_hat_ref)).*abs(indexvec).^theta);
    
        % Lawson method

        zzz_hat=zeros(size(w_hat));
        zzz_hat=u0_hat;
        tic
        for m=1:M
            M-m;
            zzzstep=zzz_hat;
            zzz_hat=kdv_onestageLawson1(zzzstep,h);
        end
        cputime_lawson(l)=toc;
        err_lawson(l)=norm(abs((zzz_hat-u_hat_ref)).*abs(indexvec).^theta);
        hl(l)=h;
    end    

    save(filename);
end