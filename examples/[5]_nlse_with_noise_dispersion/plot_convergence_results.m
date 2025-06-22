%% Here we plot the pathwise error of the numerical methods
clear all

%% Parameters
exp_settings={{'H4',2^6,2^7,1.0,2^6,2000,0.1,1}}; % Format: {'text_reg',N,Nref,mu,N_obs,Ntau}
jjj=1;
% Load parameters
text_reg=exp_settings{jjj}{1}; % Regularity parameter: options = 'Cinf', 'H2', 'H3', 'H4' 
N=exp_settings{ jjj}{2}; % Number of Fourier modes
Nref=exp_settings{jjj}{3}; % Number of Fourier modes in reference solution
mu=exp_settings{jjj}{4}; % Equation we consider is i\partial_t u= -\Delta u +mu|u|^2u
N_obs=exp_settings{jjj}{5}; % Number of observations in Wong-Zakai, i.e. \delta=T/N_obs
Ntau=exp_settings{jjj}{6}; % Number of different timesteps used to plot convergence graphs

epsilon=exp_settings{jjj}{7}; % Scaling of L^2 norm of initial data
X=exp_settings{jjj}{8}; % Number of simulation paths

T=1.0;

dataset=strcat('convergence_plot_NLSE_Tend_',num2str(T),'_',text_reg,'_N_',num2str(N),'_Nref_',num2str(Nref),'_Ntau_',num2str(Ntau),'_Nobs_',num2str(N_obs),'_epsilon_',strrep(num2str(epsilon),'.','-'),'_X_',num2str(X),'_piecewise_linear_varorig');
load(strcat('data/',dataset,'.mat'));


%% Extract sup|W(t+tau)-W(t)|
supWdiff=zeros(size(tauj));
for jj=1:max(size(tauj))
    aux=Bfunc(linspace(0,T,ceil(T/tauj(jj))+1));
    supWdiff(jj)=max(abs(aux(2:end)-aux(1:end-1)));
end

%% Convergence plot

figure(1)
clear h
loglog(tauj(:), H1_exponential(:),'^-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
hold on
%loglog(tauj(:), H1_splitting(:),'*-','linewidth',2, 'MarkerSize',10,'Color','#006400','MarkerFaceColor','white')
loglog(tauj(:), H1_relaxedCN(:),'v-','linewidth',2, 'MarkerSize',10,'Color','#FF5733','MarkerFaceColor','white')
loglog(tauj(:), H1_low_reg_new(:),'>-','linewidth',2, 'MarkerSize',10,'Color','blue','MarkerFaceColor','white')

h(1)=loglog(tauj(:), nan(size(tauj(:))),'^-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
%h(2)=loglog(tauj(:), nan(size(tauj(:))),'*-','linewidth',2, 'MarkerSize',10,'Color','#006400','MarkerFaceColor','white')
h(2)=loglog(tauj(:), nan(size(tauj(:))),'v-','linewidth',2, 'MarkerSize',10,'Color','#FF5733','MarkerFaceColor','white')
h(3)=loglog(tauj(:), nan(size(tauj(:))),'>-','linewidth',2, 'MarkerSize',10,'Color','blue','MarkerFaceColor','white')

h(4)=loglog(tauj(:),0.002*tauj(:),'--','linewidth',2,'Color','#D95319')
h(5)=loglog(tauj(:),supWdiff(:),'-.','linewidth',2,'Color','black')
loglog(tauj(:), 0.002*tauj(:).^0.5,'--','linewidth',2,'Color','#D95319')

set(gca,'FontSize',16)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 22)
if X==1
    ylabel('$\|u^{\delta,N,n}_\omega-u^{\delta,N}_\omega(T)\|_{H^1}$','Interpreter','latex', 'FontSize', 22)
else
    ylabel('$\mathbf{E}\left[\|u^{\delta,N,n}-u^{\delta,N}(T)\|_{H^1}^2\right]^{\frac{1}{2}}$','Interpreter','latex', 'FontSize', 22)
end

xlim([10^(-2.5),max(tauj(:))])
legend(h, 'Exponential','Relaxed Crank--Nicolson','SDLRI', '$O(\tau^{0.5}),O(\tau)$','$\sup_{n}|W_\omega^\delta(t_n+\tau)-W_\omega^\delta(t_n)|$','Interpreter','latex', 'FontSize', 16,'Position',[0.65 0.15 0.1 0.25])
grid on
hold off


set(gcf, 'Position',  [100, 100, 750, 600])
%exportgraphics(gcf,strcat('images/',dataset,'.pdf'),'ContentType','vector')
