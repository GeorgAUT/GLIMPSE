%% Here we test the L2 norm preservation of our symplectic method against prior work
% cubic NLSE
%
clear all
rng(10) % Fix random seed for reproducability

%% Parameters
%
text_reg='Cinf' % Regularity parameter: options = 'Cinf', 'H2', 'H3', 'H4' 
N=256; % Number of Fourier modes
Xvec=pi*(-N/2+1:N/2)'/N*2; % Spacial mesh
indexvec=(-N/2+1:N/2)'; % Vector of Fourier indices
indexvec(N/2)=1; % Correct zero index so can use to rescale Fourier coefficients
mu=1.0; % Equation we consider is i\partial_t u= -\Delta u +|u|^2u
T=20.0; % Final time

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

% Define Hamiltonian
Ham=@(u) norm(dx(u))^2+1/2*norm(conv1(u,u))^2

%% Compute the normalisation & Hamiltonian after each step
CFL=100;
tau=0.02;
M=floor(T/tau);
tau=T/M;

normalisation_first_order=zeros(1,M);
normalisation_new_symm_first_order=zeros(1,M);
normalisation_new_symm_second_order=zeros(1,M);
normalisation_brunedschratz=zeros(1,M);
normalisation_symmetric=zeros(1,M);
normalisation_strang=zeros(1,M);
normalisation_lie=zeros(1,M);
normalisation_lawson=zeros(1,M);
normalisation_res_midpoint=zeros(1,M);
hamiltonian_first_order=zeros(1,M);
hamiltonian_new_symm_first_order=zeros(1,M);
hamiltonian_new_symm_second_order=zeros(1,M);
hamiltonian_brunedschratz=zeros(1,M);
hamiltonian_symmetric=zeros(1,M);
hamiltonian_strang=zeros(1,M);
hamiltonian_lie=zeros(1,M);
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
    % Standard first order (Ostermann & Schratz)
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

    % Res-based midpoint
    wwww_hat=cubic_nls_resonance_based_symplectic(wwww_hat,mu,tau);

    normalisation_res_midpoint(j)=norm(wwww_hat);
    hamiltonian_res_midpoint(j)=Ham(wwww_hat);

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

    % Lawson
    zzz_hat=cubic_nls_symm_Lawson(zzz_hat,1,tau);

    normalisation_lawson(j)=norm(zzz_hat);
    hamiltonian_lawson(j)=Ham(zzz_hat);

    t_j(j)=j*tau; % Store time step corresponding to experiment

end

save(strcat('data/nls_structure_preservation_',text_reg,'_N_',num2str(N),'_T_',num2str(T),'_tau_',strrep(num2str(tau),'.','-'),'.mat'))

%% Plotting normalisation
figure(1)

semilogy(t_j,abs(normalisation_brunedschratz-1),'linewidth',2)
hold on
semilogy(t_j,abs(normalisation_res_midpoint-1),'linewidth',2)
semilogy(t_j,abs(normalisation_strang-1),'linewidth',2)
semilogy(t_j,abs(normalisation_lawson-1),'linewidth',2)
hold off
set(gca,'FontSize',16)
xlabel('$n\tau$','Interpreter','latex', 'FontSize', 22)
ylabel(['$\left|\|u^n\|_{L^2}-1\right|$'],'Interpreter','latex', 'FontSize', 22)

ylim([1e-16,1])
legend('Bruned \& Schratz','Res. sympl. Midpoint', 'Strang', 'Lawson','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.15 0.1 0.25])
grid on
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%% Plotting Hamiltonian
figure(2)

semilogy(t_j,abs(hamiltonian_brunedschratz-Ham(u0_hat)),'linewidth',2)
hold on
semilogy(t_j,abs(hamiltonian_res_midpoint-Ham(u0_hat)),'linewidth',2)
semilogy(t_j,abs(hamiltonian_strang-Ham(u0_hat)),'linewidth',2)
semilogy(t_j,abs(hamiltonian_lawson-Ham(u0_hat)),'linewidth',2)
hold off
set(gca,'FontSize',16)
xlabel('$n\tau$','Interpreter','latex', 'FontSize', 22)
ylabel('Absolute Error in Hamiltonian','Interpreter','latex', 'FontSize', 22)
grid on
ylim([1e-5,1])
legend('Bruned \& Schratz','Res. sympl. Midpoint', 'Strang', 'Lawson','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.15 0.1 0.25])

set(gcf, 'Position',  [100, 100, 1400, 600])