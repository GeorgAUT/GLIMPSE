%% Here we test the structure preservation of our symmetric methods against prior work
% cubic NLS equation
% The code allows for multiple sets of parameters to be specified and run
% in sequence. The data is saved and can then be plotted using
% plotting_structure_preservation.m

clear all
rng(10) % Fix random seed for reproducability

%% Operational vector
exp_settings={{'Cinf',512,40,0.02},...
    {'H2',512,40,0.02}}; % Format: {'text_reg',N,T,tau}

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
    
    u0_hat=u0_hat/norm(u0_hat); % normalisation
    
    % Define Hamiltonian
    Ham=@(u) norm(dx(u))^2+1/2*norm(conv1(u,u))^2
    %% Compute the normalisation & Hamiltonian after each step
    M=floor(T/tau);
    tau=T/M;
    
    normalisation_first_order=zeros(1,M);
    normalisation_new_symm_first_order=zeros(1,M);
    normalisation_new_symm_second_order=zeros(1,M);
    normalisation_brunedschratz=zeros(1,M);
    normalisation_symmetric=zeros(1,M);
    normalisation_strang=zeros(1,M);
    normalisation_lie=zeros(1,M);
    hamiltonian_first_order=zeros(1,M);
    hamiltonian_new_symm_first_order=zeros(1,M);
    hamiltonian_new_symm_second_order=zeros(1,M);
    hamiltonian_brunedschratz=zeros(1,M);
    hamiltonian_symmetric=zeros(1,M);
    hamiltonian_strang=zeros(1,M);
    hamiltonian_lie=zeros(1,M);
    t_j=zeros(1,M);
    
    u_hat=u0_hat;
    v_hat=u0_hat;
    w_hat=u0_hat;
    ww_hat=u0_hat;
    www_hat=u0_hat;
    z_hat=u0_hat;
    zz_hat=u0_hat;
    
    for j=1:M
        M-j
        % Standard first order
        u_hat=cubic_nls_resonance_based_first_order(u_hat,mu,tau);
        
        normalisation_first_order(j)=norm(u_hat);
        hamiltonian_first_order(j)=Ham(u_hat);
    
        % New symmetric first order
        ww_hat=cubic_nls_resonance_based_first_order_symm_new(ww_hat,mu,tau);
    
        normalisation_new_symm_first_order(j)=norm(ww_hat);
        hamiltonian_new_symm_first_order(j)=Ham(ww_hat);
    
        % New symmetric second order
        www_hat=cubic_nls_resonance_based_second_order_symm_new(www_hat,mu,tau);
    
        normalisation_new_symm_second_order(j)=norm(www_hat);
        hamiltonian_new_symm_second_order(j)=Ham(www_hat);
    
        % Symmetrised second order
        v_hat=cubic_nls_resonance_based_first_order_adjoint1(cubic_nls_resonance_based_first_order(v_hat,mu,tau/2),mu,tau/2);
        
        normalisation_symmetric(j)=norm(v_hat);
        hamiltonian_symmetric(j)=Ham(v_hat);
    
        % Standard second order
        w_hat=cubic_nls_resonance_based_second_order_BS22(w_hat,tau);
    
        normalisation_brunedschratz(j)=norm(w_hat);
        hamiltonian_brunedschratz(j)=Ham(w_hat);
    
        % Lie
        z_hat=expilaplacian(nls_nonlinear_part(z_hat,tau),tau);
    
        normalisation_lie(j)=norm(z_hat);
        hamiltonian_lie(j)=Ham(z_hat);
    
        % Strang
        zz_hat=expilaplacian(nls_nonlinear_part(expilaplacian(zz_hat,tau/2),tau),tau/2);
    
        normalisation_strang(j)=norm(zz_hat);
        hamiltonian_strang(j)=Ham(zz_hat);
    
        t_j(j)=j*tau; % Store time step corresponding to experiment
    
    end
    
    
    
    save(strcat('data/structure_preservation_',text_reg,'_N_',num2str(N),'_T_',num2str(T),'_tau_',strrep(num2str(tau),'.','-'),'.mat'))

end

