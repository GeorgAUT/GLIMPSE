%% Here we test the convergence of our symplectic method against prior work
% cubic NLSE
%
clear all
rng(10) % Fix random seed for reproducability

%% Parameters
%
text_reg='Cinf'; % Regularity parameter: options = 'Cinf', 'H2', 'H3', 'H4' 
N=2^10; % Number of Fourier modes
Xvec=pi*(-N/2+1:N/2)'/N*2; % Spacial mesh
indexvec=(-N/2+1:N/2)'; % Vector of Fourier indices
indexvec(N/2)=1; % Correct zero index so can use to rescale Fourier coefficients
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
end

u0_hat=u0_hat/norm(u0_hat); % normalisation of initial conditions
%% Reference solution computed with second order scheme from Bruned & Schratz '22
Mref=1e4;
tauref=T/Mref;

recompute=input('\n If you would like to recompute the reference sln please enter 1: ');

if recompute==1
    u_hat_ref=u0_hat;
    for m=1:Mref
        Mref-m
        
        % Standard second order
        u_hat_ref=cubic_nls_resonance_based_second_order_BS22(u_hat_ref,tauref);
    
    end
    save(strcat('data/nls_reference_solution_',text_reg,'.mat'),'u_hat_ref');
else
    load(strcat('data/nls_reference_solution_',text_reg,'.mat'),'u_hat_ref');
end

%% Compute L^2 error in time-stepping for decreasing tau
Jmax=10;

error_first_order=zeros(1,Jmax);
error_new_symm_first_order=zeros(1,Jmax);
error_new_symm_second_order=zeros(1,Jmax);
error_res_midpoint=zeros(1,Jmax);
error_brunedschratz=zeros(1,Jmax);
error_symmetric=zeros(1,Jmax);
error_strang=zeros(1,Jmax);
error_lie=zeros(1,Jmax);
error_lawson=zeros(1,Jmax);

cputime_first_order=zeros(1,Jmax);
cputime_new_symm_first_order=zeros(1,Jmax);
cputime_new_symm_second_order=zeros(1,Jmax);
cputime_res_midpoint=zeros(1,Jmax);
cputime_brunedschratz=zeros(1,Jmax);
cputime_symmetric=zeros(1,Jmax);
cputime_strang=zeros(1,Jmax);
cputime_lie=zeros(1,Jmax);
cputime_lawson=zeros(1,Jmax);
tau_j=zeros(1,Jmax);

s=1.0; % Order of Sobolev norm

parfor j=1:Jmax
    Jmax-j
    tau=1/(10*j)^1.5;
    M=floor(T/tau);
    tau=T/M;

    % Initialise each method
    u_hat=u0_hat;
    v_hat=u0_hat;
    w_hat=u0_hat;
    ww_hat=u0_hat;
    www_hat=u0_hat;
    wwww_hat=u0_hat;
    z_hat=u0_hat;
    zz_hat=u0_hat;
    zzz_hat=u0_hat;
    
    % Standard first order resonance-based scheme (Ostermann & Schratz)
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
    
    % Resonance-based midpoint rule
    tic
    for m=1:M
        wwww_hat=cubic_nls_resonance_based_symplectic(wwww_hat,1,tau);
    end
    cputime_res_midpoint(j)=toc; % CPU-time
    error_res_midpoint(j)=norm((wwww_hat-u_hat_ref).*indexvec.^s); % H^s-error
    

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

    % Lawson
    tic
    for m=1:M
        zzz_hat=cubic_nls_symm_Lawson(zzz_hat,1,tau);
    end
    cputime_lawson(j)=toc;
    error_lawson(j)=norm((zzz_hat-u_hat_ref).*indexvec.^s);

    tau_j(j)=tau; % Store time step corresponding to experiment

end    


%% Plotting
figure(1)
loglog(tau_j, error_first_order,'linewidth',2)
hold on
loglog(tau_j, error_new_symm_first_order,'linewidth',2)
loglog(tau_j, error_new_symm_second_order,'linewidth',2)
loglog(tau_j, error_brunedschratz,'linewidth',2)
loglog(tau_j, error_symmetric,'linewidth',2)
loglog(tau_j, error_res_midpoint,'linewidth',2)
loglog(tau_j, error_lie,'linewidth',2)
loglog(tau_j, error_strang,'linewidth',2)
loglog(tau_j, error_lawson,'-.','linewidth',2)
loglog(tau_j, tau_j,'--','linewidth',2,'Color','#D95319')
loglog(tau_j, tau_j.^2,'--','linewidth',2,'Color','#D95319')
set(gca,'FontSize',16)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 22)
ylabel('$H^1$-error at $T=1$','Interpreter','latex', 'FontSize', 22)
legend('Ostermann \& Schratz', 'New midpoint first order (3.12)', 'New midpoint second order (3.16)', 'Bruned \& Schratz', 'Symmetrised Res.-Based','Res. sympl. Midpoint', 'Lie', 'Strang', 'Lawson','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.15 0.1 0.25])
grid on
hold off


%% Plotting
figure(2)
loglog(tau_j, error_brunedschratz,'linewidth',2)
hold on
loglog(tau_j, error_res_midpoint,'linewidth',2)
loglog(tau_j, error_strang,'linewidth',2)
loglog(tau_j, error_lawson,'-.','linewidth',2)
loglog(tau_j, tau_j,'--','linewidth',2,'Color','#D95319')
loglog(tau_j, tau_j.^2,'--','linewidth',2,'Color','#D95319')
set(gca,'FontSize',16)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 22)
ylabel('$H^1$-error at $T=1$','Interpreter','latex', 'FontSize', 22)
legend('Bruned \& Schratz','Res. sympl. Midpoint', 'Strang', 'Lawson','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.15 0.1 0.25])
grid on
hold off


set(gcf, 'Position',  [100, 100, 700, 600])
%exportgraphics(gcf,'images/convergence_plot_NLSE_H1_Cinf.pdf','ContentType','vector')

%% Time plots Plotting
figure(3)
loglog(cputime_brunedschratz,error_brunedschratz ,'linewidth',2)
hold on
loglog(cputime_res_midpoint,error_res_midpoint,'linewidth',2)
loglog(cputime_strang,error_strang,'linewidth',2)
loglog(cputime_lawson,error_lawson,'-.','linewidth',2)
%loglog(tau_j, tau_j,'--','linewidth',2,'Color','#D95319')
%loglog(tau_j, tau_j.^2,'--','linewidth',2,'Color','#D95319')
set(gca,'FontSize',16)
xlabel('CPU time (sec)','Interpreter','latex', 'FontSize', 22)
ylabel('$H^1$-error at $T=1$','Interpreter','latex', 'FontSize', 22)
legend('Bruned \& Schratz','Res. sympl. Midpoint', 'Strang', 'Lawson','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.15 0.1 0.25])
grid on
hold off


set(gcf, 'Position',  [100, 100, 700, 600])
%exportgraphics(gcf,'images/cpu_time_plot_NLSE_H1_Cinf.pdf','ContentType','vector')



