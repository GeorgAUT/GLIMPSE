%% Here we test the structure preservation of our symmetric methods against prior work
% cubic NLS equation in 2D
% The code allows for multiple sets of parameters to be specified and run
% in sequence. The data is saved and can then be plotted using
% plotting_structure_preservation.m

clear all
rng(10) % Fix random seed for reproducability
norm2d=@(M) sqrt(sum(sum(abs(M).^2))); % Define l2 norm for matrix valued entries

%% Operational vector
exp_settings={{'Cinf',16,100,0.02,1},...
    {'H2',16,100,0.02,-1}}; % Format: {'text_reg',N,T,tau,mu}

for jjj=1:size(exp_settings,2)

    text_reg=exp_settings{jjj}{1};
    N=exp_settings{jjj}{2};
    T=exp_settings{jjj}{3};
    tau=exp_settings{jjj}{4};
    mu=exp_settings{jjj}{5};% Equation we consider is i\partial_t u= -\Delta u +|u|^2u
    %% Parameters
    %
    xvec=pi*(-N/2+1:N/2)'/N*2; % Spacial mesh
    yvec=pi*(-N/2+1:N/2)'/N*2; % Spacial mesh
    
    [Xvec,Yvec]=meshgrid(xvec,yvec);
    laplace_index=(((-N/2+1:N/2).^2)'+(-N/2+1:N/2).^2); % Vector of Fourier indices
    jap_brack=(1+laplace_index).^(1/2); % Japanese bracket for Sobolev norms
    
    
    %% Initial conditions
    %
    switch text_reg
        case 'Cinf'
            u0_hat=fft2pi(cos(Yvec)./(2+sin(Yvec)).*sin(Xvec)./(2+cos(Xvec))); %smooth initial data
        case 'H2'
            theta=2.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./jap_brack.^theta;
        case 'H3'
            theta=3.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./jap_brackc.^theta;
        case 'H4'
            theta=4.0; %regularity parameter for initial condition
            u0_hat=rand(size(Xvec))+i*rand(size(Xvec));
            u0_hat=u0_hat./jap_brack.^theta;
    end
    
    u0_hat=u0_hat/norm2d(u0_hat); % normalisation
    
    % Define Hamiltonian
    Ham=@(u) norm2d(dx2d(u))^2+norm2d(dy2d(u))^2+mu/2*norm2d(conv2d(u,u))^2;
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
    H2norm_first_order=zeros(1,M);
    H2norm_new_symm_first_order=zeros(1,M);
    H2norm_new_symm_second_order=zeros(1,M);
    H2norm_brunedschratz=zeros(1,M);
    H2norm_symmetric=zeros(1,M);
    H2norm_strang=zeros(1,M);
    H2norm_lie=zeros(1,M);
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
        u_hat=cubic_nls_2d_resonance_based_first_order(u_hat,mu,tau);
        
        normalisation_first_order(j)=norm2d(u_hat);
        H2norm_first_order(j)=norm2d(u_hat.*jap_brack.^2);
        hamiltonian_first_order(j)=Ham(u_hat);
    
        % New symmetric first order
        ww_hat=cubic_nls_2d_resonance_based_first_order_symm_new(ww_hat,mu,tau);
    
        normalisation_new_symm_first_order(j)=norm2d(ww_hat);
        H2norm_new_symm_first_order(j)=norm2d(ww_hat.*jap_brack.^2);
        hamiltonian_new_symm_first_order(j)=Ham(ww_hat);
    
        % New symmetric second order
        www_hat=cubic_nls_2d_resonance_based_second_order_symm_new(www_hat,mu,tau);
    
        normalisation_new_symm_second_order(j)=norm2d(www_hat);
        H2norm_new_symm_second_order(j)=norm2d(www_hat.*jap_brack.^2);
        hamiltonian_new_symm_second_order(j)=Ham(www_hat);
    
        % Symmetrised second order
        v_hat=cubic_nls_2d_resonance_based_first_order_adjoint(cubic_nls_2d_resonance_based_first_order(v_hat,mu,tau/2),mu,tau/2);
        
        normalisation_symmetric(j)=norm2d(v_hat);
        H2norm_symmetric(j)=norm2d(v_hat.*jap_brack.^2);
        hamiltonian_symmetric(j)=Ham(v_hat);
    
        % Standard second order
        w_hat=cubic_nls_2d_resonance_based_second_order_BS22(w_hat,mu,tau);
    
        normalisation_brunedschratz(j)=norm2d(w_hat);
        H2norm_brunedschratz(j)=norm2d(w_hat.*jap_brack.^2);
        hamiltonian_brunedschratz(j)=Ham(w_hat);
    
        % Lie
        z_hat=expilaplacian2d(nls_nonlinear_part2d(z_hat,tau),tau);
    
        normalisation_lie(j)=norm2d(z_hat);
        H2norm_lie(j)=norm2d(z_hat.*jap_brack.^2);
        hamiltonian_lie(j)=Ham(z_hat);
    
        % Strang
        zz_hat=expilaplacian2d(nls_nonlinear_part2d(expilaplacian2d(zz_hat,tau/2),tau),tau/2);
    
        normalisation_strang(j)=norm2d(zz_hat);
        H2norm_strang(j)=norm2d(zz_hat.*jap_brack.^2);
        hamiltonian_strang(j)=Ham(zz_hat);
    
        t_j(j)=j*tau; % Store time step corresponding to experiment

        if mod(j,100)==0 % Intermediate saving
            save(strcat('data/structure_preservation_2d_',text_reg,'_N_',num2str(N),'_T_',num2str(T),'_tau_',strrep(num2str(tau),'.','-'),'_mu_',strrep(num2str(mu),'.','-'),'.mat'))
        end
        
    end
    
    
    
    save(strcat('data/structure_preservation_2d_',text_reg,'_N_',num2str(N),'_T_',num2str(T),'_tau_',strrep(num2str(tau),'.','-'),'_mu_',strrep(num2str(mu),'.','-'),'.mat'))

end

