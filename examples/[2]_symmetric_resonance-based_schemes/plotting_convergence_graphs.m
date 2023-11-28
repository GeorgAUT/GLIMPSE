%% Plotting of convergence graphs
% Data can be produced with convergence_experiments_nls_multiple_runs.m

clear all
rng(10) % Fix random seed for reproducability

%% Plotting

dataset='H2_Nref_256_N_64_T_1_Jmax_10';
clear h
load(strcat('data/convergence_data_nlse_',dataset,'.mat'))
figure(3)

loglog(tau_j, error_first_order,'linewidth',2,'Color','#0072bd','MarkerFaceColor','white')
hold on
loglog(tau_j, error_new_symm_first_order,'linewidth',2,'Color','#7e2f8e','MarkerFaceColor','white')
loglog(tau_j, error_new_symm_second_order,'linewidth',2,'Color','#edb120','MarkerFaceColor','white')
loglog(tau_j, error_brunedschratz,'linewidth',2,'Color','black','MarkerFaceColor','white')
loglog(tau_j, error_symmetric,'linewidth',2,'Color','#c1930d','MarkerFaceColor','white')
loglog(tau_j, error_lie,'linewidth',2,'Color','#ee73c4','MarkerFaceColor','white')
loglog(tau_j, error_strang,'linewidth',2,'Color','#1d7d01','MarkerFaceColor','white')

step=10;
loglog(tau_j(1:step:end), error_first_order(1:step:end),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
loglog(tau_j(1:step:end), error_new_symm_first_order(1:step:end),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
loglog(tau_j(1:step:end), error_new_symm_second_order(1:step:end),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
loglog(tau_j(1:step:end), error_brunedschratz(1:step:end),'d','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
loglog(tau_j(1:step:end), error_symmetric(1:step:end),'^','linewidth',2, 'MarkerSize',10,'Color','#c1930d','MarkerFaceColor','white')
loglog(tau_j(1:step:end), error_lie(1:step:end),'>','linewidth',2, 'MarkerSize',10,'Color','#ee73c4','MarkerFaceColor','white')
loglog(tau_j(1:step:end), error_strang(1:step:end),'square','linewidth',2, 'MarkerSize',10,'Color','#1d7d01','MarkerFaceColor','white')


h(1)=loglog(tau_j, nan(size(tau_j)),'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(2)=loglog(tau_j, nan(size(tau_j)),'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(3)=loglog(tau_j, nan(size(tau_j)),'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
h(4)=loglog(tau_j, nan(size(tau_j)),'d-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(5)=loglog(tau_j, nan(size(tau_j)),'^-','linewidth',2, 'MarkerSize',10,'Color','#c1930d','MarkerFaceColor','white')
h(6)=loglog(tau_j, nan(size(tau_j)),'>-','linewidth',2, 'MarkerSize',10,'Color','#ee73c4','MarkerFaceColor','white')
h(7)=loglog(tau_j, nan(size(tau_j)),'square-','linewidth',2, 'MarkerSize',10,'Color','#1d7d01','MarkerFaceColor','white')

h(8)=loglog(tau_j, tau_j,'--','linewidth',2,'Color','#D95319')
loglog(tau_j, tau_j.^2,'--','linewidth',2,'Color','#D95319')
set(gca,'FontSize',16)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 22)
ylabel('$H^1$-error at $T=1$','Interpreter','latex', 'FontSize', 22)

xlim([min(tau_j),max(tau_j)])
%xlim([1e-3,1])
legend(h,'Ostermann \& Schratz', 'New midpoint first order (3.14)', 'New midpoint second order (3.15)', 'Bruned \& Schratz', 'Alama Bronsard', 'Lie', 'Strang','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 16,'Position',[0.65 0.15 0.1 0.25])
grid on
hold off


set(gcf, 'Position',  [100, 100, 750, 600])
%exportgraphics(gcf,strcat('images/convergence_plot_NLSE_H1_',dataset,'.pdf'),'ContentType','vector')
