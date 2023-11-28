%% Plotting of structure preservation results
% Data can be generated with structure_preservation_multiple_runs.m

clear all
addpath('modules/')

%% Plot for examples section (energy)
figure(1)
clear h;
nstep=2500;
dataset='Cinf_N_512_T_40_tau_0-02'
load(strcat('data/structure_preservation_',dataset,'.mat'))
semilogy(t_j,ones(size(t_j)),'-.','linewidth',1.5, 'Color','red')
hold on
semilogy(t_j,abs(hamiltonian_first_order-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')
semilogy(t_j(1:nstep:end),abs(hamiltonian_first_order(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')

semilogy(t_j,abs(hamiltonian_new_symm_first_order-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')%'#8034eb')
semilogy(t_j(1:nstep:end),abs(hamiltonian_new_symm_first_order(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')

semilogy(t_j,abs(hamiltonian_new_symm_second_order-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')
semilogy(t_j(1:nstep:end),abs(hamiltonian_new_symm_second_order(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')

semilogy(t_j,abs(hamiltonian_brunedschratz-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','black')
semilogy(t_j(1:nstep:end),abs(hamiltonian_brunedschratz(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'d','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')

semilogy(t_j,abs(hamiltonian_symmetric-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#c1930d')
semilogy(t_j(1:nstep:end),abs(hamiltonian_symmetric(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'^','linewidth',2, 'MarkerSize',10,'Color','#c1930d','MarkerFaceColor','white')

semilogy(t_j,abs(hamiltonian_lie-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#ee73c4')
semilogy(t_j(1:nstep:end),abs(hamiltonian_lie(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'>','linewidth',2, 'MarkerSize',10,'Color','#ee73c4','MarkerFaceColor','white')

semilogy(t_j,abs(hamiltonian_strang-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#1d7d01')
semilogy(t_j(1:nstep:end),abs(hamiltonian_strang(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'square','linewidth',2, 'MarkerSize',10,'Color','#1d7d01','MarkerFaceColor','white')

h(1)=semilogy(t_j, nan(size(t_j)),'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(2)=semilogy(t_j, nan(size(t_j)),'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(3)=semilogy(t_j, nan(size(t_j)),'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
h(4)=semilogy(t_j, nan(size(t_j)),'d-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(5)=semilogy(t_j, nan(size(t_j)),'^-','linewidth',2, 'MarkerSize',10,'Color','#c1930d','MarkerFaceColor','white')
h(6)=semilogy(t_j, nan(size(t_j)),'>-','linewidth',2, 'MarkerSize',10,'Color','#ee73c4','MarkerFaceColor','white')
h(7)=semilogy(t_j, nan(size(t_j)),'square-','linewidth',2, 'MarkerSize',10,'Color','#1d7d01','MarkerFaceColor','white')

hold off

set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Rel. error in energy','Interpreter','latex', 'FontSize', 22)

legend(h,'Ostermann \& Schratz', 'New midpoint first order (3.14)', 'New midpoint second order (3.15)', 'Bruned \& Schratz', 'Alama Bronsard', 'Lie', 'Strang','Interpreter','latex', 'FontSize', 16,'Position',[0.74 0.68 0.1 0.2])
grid on


ylim([1e-10,1e2])
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/energy_graph_nlse_',dataset,'.pdf'),'ContentType','vector')

%% Plot for examples section (normalisation)
figure(3)
clear h;
nstep=2500;
load(strcat('data/structure_preservation_',dataset,'.mat'))
semilogy(t_j,ones(size(t_j)),'-.','linewidth',1.5, 'Color','red')
hold on
semilogy(t_j,abs(normalisation_first_order-1),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')
semilogy(t_j(1:nstep:end),abs(normalisation_first_order(1:nstep:end)-1),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')

semilogy(t_j,abs(normalisation_new_symm_first_order-1),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')%'#8034eb')
semilogy(t_j(1:nstep:end),abs(normalisation_new_symm_first_order(1:nstep:end)-1),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')

semilogy(t_j,abs(normalisation_new_symm_second_order-1),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')
semilogy(t_j(1:nstep:end),abs(normalisation_new_symm_second_order(1:nstep:end)-1),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')

semilogy(t_j,abs(normalisation_brunedschratz-1),'-','linewidth',2, 'MarkerSize',10,'Color','black')
semilogy(t_j(1:nstep:end),abs(normalisation_brunedschratz(1:nstep:end)-1),'d','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')

semilogy(t_j,abs(normalisation_symmetric-1),'-','linewidth',2, 'MarkerSize',10,'Color','#c1930d')
semilogy(t_j(1:nstep:end),abs(normalisation_symmetric(1:nstep:end)-1),'^','linewidth',2, 'MarkerSize',10,'Color','#c1930d','MarkerFaceColor','white')

semilogy(t_j,abs(normalisation_lie-1),'-','linewidth',2, 'MarkerSize',10,'Color','#ee73c4')
semilogy(t_j(1:nstep:end),abs(normalisation_lie(1:nstep:end)-1),'>','linewidth',2, 'MarkerSize',10,'Color','#ee73c4','MarkerFaceColor','white')

semilogy(t_j,abs(normalisation_strang-1),'-','linewidth',2, 'MarkerSize',10,'Color','#1d7d01')
semilogy(t_j(1:nstep:end),abs(normalisation_strang(1:nstep:end)-1),'square','linewidth',2, 'MarkerSize',10,'Color','#1d7d01','MarkerFaceColor','white')

h(1)=semilogy(t_j, nan(size(t_j)),'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(2)=semilogy(t_j, nan(size(t_j)),'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(3)=semilogy(t_j, nan(size(t_j)),'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
h(4)=semilogy(t_j, nan(size(t_j)),'d-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')
h(5)=semilogy(t_j, nan(size(t_j)),'^-','linewidth',2, 'MarkerSize',10,'Color','#c1930d','MarkerFaceColor','white')
h(6)=semilogy(t_j, nan(size(t_j)),'>-','linewidth',2, 'MarkerSize',10,'Color','#ee73c4','MarkerFaceColor','white')
h(7)=semilogy(t_j, nan(size(t_j)),'square-','linewidth',2, 'MarkerSize',10,'Color','#1d7d01','MarkerFaceColor','white')


hold off

set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Rel. error in $\|u\|_{L^2}^2$','Interpreter','latex', 'FontSize', 22)

legend(h,'Ostermann \& Schratz', 'New midpoint first order (3.14)', 'New midpoint second order (3.15)', 'Bruned \& Schratz', 'Alama Bronsard', 'Lie', 'Strang','Interpreter','latex', 'FontSize', 16,'Position',[0.74 0.68 0.1 0.2])
grid on


ylim([1e-16,1e2])
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/normalisaiton_graph_nlse_',dataset,'.pdf'),'ContentType','vector')
