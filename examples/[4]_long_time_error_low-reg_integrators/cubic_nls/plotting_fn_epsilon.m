%% Plotting the results of epsilon dependence in our new schemes
% Data can be produced with long_time_fn_epsilon_dependence.m
%
clear h
load('data/long_time_function_of_epsilon_tau_0_05_H2_data_Ntest_1024_Nref_4096_v1.mat')
figure(1)
h(1)=loglog(epsilon_j,error_res,'v-','linewidth',3, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
h(2)=loglog(epsilon_j,error_res_long_time,'|-','linewidth',3, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
h(3)=loglog(epsilon_j,error_res_long_time_second_order,'<-','linewidth',3, 'MarkerSize',10,'Color','#c334eb','MarkerFaceColor','white')
h(4)=loglog(epsilon_j,error_lie,'o-','linewidth',3, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(5)=loglog(epsilon_j,error_strang,'^-','linewidth',3, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')

loglog(epsilon_j,epsilon_j.^2*2*1e-1,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')
h(6)=loglog(nan,nan,'--','linewidth',2.0, 'MarkerSize',10,'Color','#D95319')


set(gca,'FontSize',20)
xlabel('$\varepsilon$','Interpreter','latex', 'FontSize', 26)
ylabel('$H^1$-error at $T=1/\varepsilon^2$','Interpreter','latex', 'FontSize', 26)
legend(h,"Ostermann \& Schratz '18",'NRLI1','NRSLI2','Lie','Strang','$O(\varepsilon^2)$','Interpreter','latex', 'FontSize', 20,'Position',[0.66 0.14 0.1 0.25])
grid on
set(gcf, 'Position',  [100, 100, 500, 400])

