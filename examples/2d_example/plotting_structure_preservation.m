clear all

%% Select the dataset to be plotted

dataset=input("\n\n Which dataset would you like to plot, e.g. 'Cinf_N_16_T_100_tau_0-02_mu_1' (please don't forget the '' marks):\n");

%% Load dataset

load(strcat('data/structure_preservation_2d_',dataset,'.mat'))
indexvec=jap_brack;

%% Plot H2 norm
figure(1)
clear h;
nstep=2500;
semilogy(t_j,ones(size(t_j)),'-.','linewidth',1.5, 'Color','red')
hold on
semilogy(t_j,abs(H2norm_brunedschratz-norm(u0_hat.*indexvec.^2))/norm(u0_hat.*indexvec.^2),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')

semilogy(t_j(1:nstep:end),abs(H2norm_brunedschratz(1:nstep:end)-norm(u0_hat.*indexvec.^2))/norm(u0_hat.*indexvec.^2),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
semilogy(t_j,abs(H2norm_strang-norm(u0_hat.*indexvec.^2))/norm(u0_hat.*indexvec.^2),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')%'#8034eb')
semilogy(t_j(1:nstep:end),abs(H2norm_strang(1:nstep:end)-norm(u0_hat.*indexvec.^2))/norm(u0_hat.*indexvec.^2),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
semilogy(t_j,abs(H2norm_new_symm_second_order-norm(u0_hat.*indexvec.^2))/norm(u0_hat.*indexvec.^2),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')

semilogy(t_j(1:nstep:end),abs(H2norm_new_symm_second_order(1:nstep:end)-norm(u0_hat.*indexvec.^2))/norm(u0_hat.*indexvec.^2),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')

h(1)=semilogy(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(2)=semilogy(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(3)=semilogy(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')

hold off

set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Rel. error in $H^2$ norm','Interpreter','latex', 'FontSize', 22)

legend(h,'New symmetric method','Bruned \& Schratz 2022','Strang','Interpreter','latex', 'FontSize', 16,'Position',[0.77 0.15 0.1 0.2])
grid on


ylim([1e-10,1e2])
hold off

% Store the graphic as pdf vector
set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/H2_graph_',dimensioninput,'_',dataset,'.pdf'),'ContentType','vector')

%% Plot normalisation
figure(2)
clear h;
nstep=2500;
semilogy(t_j,abs(normalisation_brunedschratz-1),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')
hold on
semilogy(t_j(1:nstep:end),abs(normalisation_brunedschratz(1:nstep:end)-1),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
semilogy(t_j,abs(normalisation_strang-1),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')%'#8034eb')
semilogy(t_j(1:nstep:end),abs(normalisation_strang(1:nstep:end)-1),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
semilogy(t_j,abs(normalisation_new_symm_second_order-1),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')

semilogy(t_j(1:nstep:end),abs(normalisation_new_symm_second_order(1:nstep:end)-1),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')

h(1)=semilogy(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(2)=semilogy(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(3)=semilogy(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')

hold off

tau
tau*N^2
set(gca,'FontSize',16)
%title('Strang, N=64')
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Error in normalisation','Interpreter','latex', 'FontSize', 22)

legend(h,'New symmetric method','Bruned \& Schratz 2022','Strang','Interpreter','latex', 'FontSize', 16,'Position',[0.77 0.68 0.1 0.2])
grid on


ylim([1e-16,1e1])
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/normalisation_graph_intro_',dimensioninput,'_',dataset,'.pdf'),'ContentType','vector')


%% Plot Hamiltonian
figure(3)
clear h;
nstep=2500;
semilogy(t_j,ones(size(t_j)),'-.','linewidth',1.5, 'Color','red')
hold on
semilogy(t_j,abs(hamiltonian_brunedschratz-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')

semilogy(t_j(1:nstep:end),abs(hamiltonian_brunedschratz(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
semilogy(t_j,abs(hamiltonian_strang-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')%'#8034eb')
semilogy(t_j(1:nstep:end),abs(hamiltonian_strang(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
semilogy(t_j,abs(hamiltonian_new_symm_second_order-Ham(u0_hat))/Ham(u0_hat),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')

semilogy(t_j(1:nstep:end),abs(hamiltonian_new_symm_second_order(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')

h(1)=semilogy(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(2)=semilogy(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(3)=semilogy(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')

hold off

set(gca,'FontSize',16)
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Rel. error in energy','Interpreter','latex', 'FontSize', 22)

legend(h,'New symmetric method','Bruned \& Schratz 2022','Strang','Interpreter','latex', 'FontSize', 16,'Position',[0.77 0.15 0.1 0.2])
grid on


ylim([1e-10,1e2])
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/energy_graph_',dimensioninput,'_',dataset,'.pdf'),'ContentType','vector')

