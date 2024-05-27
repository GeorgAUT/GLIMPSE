%% Plotting the results
load('data/error_evolution_eps_0-1_H2_data_Ntest_8192_Nref_16384.mat')
clear h

figure(1)
h(1)=semilogy(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
hold on
h(2)=semilogy(nan,nan,'|-','linewidth',2, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
h(3)=semilogy(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#c334eb','MarkerFaceColor','white')
h(4)=semilogy(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(5)=semilogy(nan,nan,'^-','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')

nstep=1;

semilogy(t_vec(1:nstep:end),error_res(1:nstep:end),'-','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_res_long_time,'-','linewidth',2, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_res_long_time_second_order,'-','linewidth',2, 'MarkerSize',10,'Color','#c334eb','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_lie(1:nstep:end),'-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_strang(1:nstep:end),'-','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')

nstep=100;

semilogy(t_vec(1:nstep:end),error_res(1:nstep:end),'v','linewidth',2, 'MarkerSize',10,'Color','#77AC30','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_res_long_time(1:nstep:end),'|','linewidth',2, 'MarkerSize',10,'Color','#8034eb','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_res_long_time_second_order(1:nstep:end),'<','linewidth',2, 'MarkerSize',10,'Color','#c334eb','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_lie(1:nstep:end),'o','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
semilogy(t_vec(1:nstep:end),error_strang(1:nstep:end),'^','linewidth',2, 'MarkerSize',10,'Color','#3498eb','MarkerFaceColor','white')

ylim([1e-5,0.25])
set(gca,'FontSize',20)
xlabel('$t$','Interpreter','latex', 'FontSize', 24)
ylabel('$H^1$-error at $t$','Interpreter','latex', 'FontSize', 24)
legend(h,"Ostermann \& Schratz '18",'NRLI1','NRSLI2','Lie','Strang','Interpreter','latex', 'FontSize', 18,'Position',[0.66 0.14 0.1 0.25])
grid on
set(gcf, 'Position',  [100, 100, 750, 600])