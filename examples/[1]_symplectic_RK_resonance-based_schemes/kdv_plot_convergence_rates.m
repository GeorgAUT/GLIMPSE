%% Plotting of convergence results
% Data can be generated with kdv_convergence_experiments_multiple_runs.m
%
clear all
addpath('modules/')
rng(10) % Fix random seed for reproducability

%% Plotting

dataset='H4_Nref_1024_N_64_T_1_Jmax_20';
clear h

load(strcat('data/kdv_convergence_data_',dataset,'.mat'))
figure(5)

err_brunedschratz=err_brunedschratz/10;
err_strang=err_strang/10;
err_lawson=err_lawson/10;
err_res_midpoint=err_res_midpoint/10;

step=20;
tend=100%500;
loglog(hl, err_brunedschratz,'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
hold on
loglog(hl, err_strang,'-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
loglog(hl, err_lawson,'-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
loglog(hl, err_res_midpoint,'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')

loglog(hl(1:step:end), err_brunedschratz(1:step:end),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
loglog(hl(1:step:end), err_strang(1:step:end),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
loglog(hl(1:step:end), err_lawson(1:step:end),'d','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
loglog(hl(1:step:end), err_res_midpoint(1:step:end),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')

h(1)=loglog(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(2)=loglog(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(3)=loglog(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
h(4)=loglog(nan,nan,'d-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')

h(5)=loglog(hl, hl*30,'--','linewidth',2,'Color','#D95319')
loglog(hl, hl.^2*30,'--','linewidth',2,'Color','#D95319')
set(gca,'FontSize',16)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 22)
ylabel('$H^1$-error at $T=1$','Interpreter','latex', 'FontSize', 22)
legend(h,'Resonance-based midpoint','Bruned \& Schratz', 'Strang', 'Lawson','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 16,'Position',[0.68 0.14 0.1 0.25])
grid on
hold off
xlim([min(hl),max(hl)])
ylim([1e-6,1e1])
set(gcf, 'Position',  [100, 100, 700, 600])
%exportgraphics(gcf,strcat('images/convergence_plot_kdv_scaled_H1_',dataset,'.pdf'),'ContentType','vector')

