%% Plotting the results
clear h
load('data/long_time_function_of_epsilon_tau_0_25_H1_data_Ntest_256_Nref_512_v2.mat')
figure(3)
h(1)=loglog(epsilon_j,error_res,'v-','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
h(2)=loglog(epsilon_j,error_sym_res,'d-','linewidth',2, 'MarkerSize',10,'Color','#0072BD','MarkerFaceColor','white')
h(3)=loglog(epsilon_j,error_lie,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(4)=loglog(epsilon_j,error_strang,'^-','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')

loglog(epsilon_j,epsilon_j*50*1e-1,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')
h(5)=loglog(nan,nan,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')


set(gca,'FontSize',20)
xlabel('$\varepsilon$','Interpreter','latex', 'FontSize', 24)
ylabel('$H^1$-error at $T=1/\varepsilon$','Interpreter','latex', 'FontSize', 24)
legend(h,'LI1','SLI2','Lie','Strang','$O(\varepsilon)$','Interpreter','latex', 'FontSize', 18,'Position',[0.77 0.13 0.1 0.25])
grid on
xlim([min(epsilon_j),1])
set(gcf, 'Position',  [100, 100, 750, 600])