%% Plotting the results
clear h
load('data/long_time_function_of_tau_epsilon_0-5_H2_data_N_256_v2')
figure(1)
h(1)=loglog(tau_j(2:end),error_res(2:end),'v-','linewidth',3, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
h(2)=loglog(tau_j(5:end),error_res_long_time(5:end),'|-','linewidth',3, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
h(3)=loglog(tau_j,error_res_long_time_second_order,'<-','linewidth',3, 'MarkerSize',10,'Color','#c334eb','MarkerFaceColor','white')
h(4)=loglog(tau_j,error_lie,'o-','linewidth',3, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(5)=loglog(tau_j,error_strang,'^-','linewidth',3, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')
loglog(tau_j,error_res_long_time_second_order,'<-','linewidth',3, 'MarkerSize',10,'Color','#c334eb','MarkerFaceColor','white')
loglog(tau_j(5:end),error_res_long_time(5:end),'|-','linewidth',3, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')

%loglog(tau_j,error_sym_res,'d-','linewidth',2, 'MarkerSize',10,'Color','#0072BD','MarkerFaceColor','white')
loglog(tau_j,tau_j*1e+2,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')
loglog(tau_j,tau_j.^2*.3*1e+1,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')
h(6)=loglog(nan,nan,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')

set(gca,'FontSize',20)
xlabel('$\tau$','Interpreter','latex', 'FontSize', 26)
ylabel('$H^1$-error at $T=1/\varepsilon^2$','Interpreter','latex', 'FontSize', 26)
legend(h,"Ostermann \& Schratz '18",'NRLI1','NRSLI2','Lie','Strang','$O(\tau),O(\tau^2)$','Interpreter','latex', 'FontSize', 20,'Position',[0.665 0.135 0.1 0.25])
grid on
%xlim([0.75*1e-2,1e-1])
set(gcf, 'Position',  [100, 100, 500, 400])