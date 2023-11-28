%% Plotting the results
clear h
N=256;
load(strcat('data/long_time_function_of_tau_epsilon_0-1_H1_data_N_',num2str(N),'.mat'))
figure(1)
h(1)=loglog(tau_j,error_res,'v-','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
h(2)=loglog(tau_j,error_sym_res,'d-','linewidth',2, 'MarkerSize',10,'Color','#0072BD','MarkerFaceColor','white')
h(3)=loglog(tau_j,error_lie,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(4)=loglog(tau_j,error_strang,'^-','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')

loglog(tau_j,tau_j*1e+0,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')
loglog(tau_j,tau_j.^2*0.7*1e-1,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')
h(5)=loglog(nan,nan,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')

set(gca,'FontSize',20)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 24)
ylabel('$H^1$-error at $T=1/\varepsilon$','Interpreter','latex', 'FontSize', 24)
legend(h,'LI1','SLI2','Lie','Strang','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 18,'Position',[0.72 0.15 0.1 0.25])
grid on
set(gcf, 'Position',  [100, 100, 750, 600])
