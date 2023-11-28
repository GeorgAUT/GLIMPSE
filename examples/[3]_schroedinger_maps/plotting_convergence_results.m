%% Plotting of the convergence behaviour of our method in comparison to prior work
% To produce data please use comparison_SM_methods.m
%
load('data/convergence_results_rough_T-1_N_256_v4.mat')
%% Plot the results (convergence rate)
figure(1)
graph(1)=loglog(hll(1:end),err_be(1:end),'v-','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
graph(2)=loglog(hll(1:end),err_chag(1:end),'d-','linewidth',2, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
graph(3)=loglog(hll(1:end),err_interpolatory_magnus_O2(1:end),'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
graph(4)=loglog(hll(1:end),err_interpolatory_magnus_O4(1:end),'^-','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')
graph(5)=loglog(hll(1:end),err_low_reg(1:end),'x-','linewidth',2, 'MarkerSize',10,'Color','#ebb734','MarkerFaceColor','white')
loglog(hll(1:end),hll(1:end).^0.5*1e-1*0.95,'--','linewidth',1.0, 'MarkerSize',10,'Color','#D95319')
graph(6)=loglog(nan,nan,'--','linewidth',1.0, 'MarkerSize',10,'Color','#D95319')

set(gca,'FontSize',16)
%xlim([1e-3,1])
%ylim([1e-3,2*1e0])
xlabel('$h$','Interpreter','latex', 'FontSize', 22)
ylabel('$L^2$-error at $T=1$','Interpreter','latex', 'FontSize', 22)
legend(graph,'Xie et al. 2020 - order 1','Xie et al. 2020 - order 2','Scheme A - order 2','Scheme A - order 4','Scheme B','$O(h^{1/2})$','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.15 0.1 0.25])
grid on
set(gcf, 'Position',  [100, 100, 750, 600])

%% Plot the results (convergence time)

clear graph
figure(2)
graph(1)=loglog(err_be(1:end),time_be(1:end),'v-','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
graph(2)=loglog(err_chag(1:end),time_chag(1:end),'d-','linewidth',2, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
graph(3)=loglog(err_interpolatory_magnus_O2(1:end),time_interpolatory_magnus_O2(1:end),'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
graph(4)=loglog(err_interpolatory_magnus_O4(1:end),time_interpolatory_magnus_O4(1:end),'^-','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')
graph(5)=loglog(err_low_reg(1:end),time_low_reg(1:end),'x-','linewidth',2, 'MarkerSize',10,'Color','#ebb734','MarkerFaceColor','white')

set(gca,'FontSize',16)
ylabel('CPU-time (sec)','Interpreter','latex', 'FontSize', 22)
xlabel('$L^2$-error at $T=1$','Interpreter','latex', 'FontSize', 22)
legend(graph,'Xie et al. 2020 - order 1','Xie et al. 2020 - order 2','Scheme A - order 2','Scheme A - order 4','Scheme B','Interpreter','latex', 'FontSize', 16,'Position',[0.70 0.66 0.1 0.25])
grid on
set(gcf, 'Position',  [100, 100, 750, 600])