%% Plotting of structure preservation results
% Data can be generated with kdv_structure_preservation_multiple_runs.m
%
clear all

dataset='scaled_Cinf_N_64_T_40_tau_0-02'
load(strcat('data/kdv_structure_preservation_',dataset,'.mat'))
%% Plot Hamiltonian
figure(1)
clear h;
nstep=2500;
semilogy(t_j,ones(size(t_j)),'-.','linewidth',1.5, 'Color','red')
hold on
semilogy(t_j,abs((hamiltonian_brunedschratz-Ham(u0_hat))/Ham(u0_hat)),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')
semilogy(t_j(1:nstep:end),abs((hamiltonian_brunedschratz(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat)),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
semilogy(t_j,abs((hamiltonian_strang-Ham(u0_hat))/Ham(u0_hat)),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')%'#8034eb')
semilogy(t_j(1:nstep:end),abs((hamiltonian_strang(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat)),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
semilogy(t_j,abs((hamiltonian_lawson-Ham(u0_hat))/Ham(u0_hat)),'-','linewidth',2, 'MarkerSize',10,'Color','black')%'#8034eb')
semilogy(t_j(1:nstep:end),abs((hamiltonian_lawson(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat)),'d','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')

semilogy(t_j,abs((hamiltonian_res_midpoint-Ham(u0_hat))/Ham(u0_hat)),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')
semilogy(t_j(1:nstep:end),abs((hamiltonian_res_midpoint(1:nstep:end)-Ham(u0_hat))/Ham(u0_hat)),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')


h(1)=semilogy(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(2)=semilogy(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(3)=semilogy(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
h(4)=semilogy(nan,nan,'d-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')

hold off

tau
tau*N^2
set(gca,'FontSize',16)
%title('Strang, N=64')
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Rel. error in $\mathcal{H}^{[KdV]}(u)$','Interpreter','latex', 'FontSize', 22)

legend(h,'Resonance-based midpoint','Bruned \& Schratz', 'Strang', 'Lawson','Interpreter','latex', 'FontSize', 16,'Position',[0.77 0.5 0.1 0.2])
grid on


ylim([1e-10,1e2])
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/kdv_energy_preservation_',dataset,'.pdf'),'ContentType','vector')

%% Plot momentum
figure(2)
clear h;
nstep=2500;
momentum_ref=(sum(ifftpi(u0_hat).^2))/N;
semilogy(t_j,ones(size(t_j)),'-.','linewidth',1.5, 'Color','red')
hold on
semilogy(t_j,abs((momentum_res_midpoint-momentum_ref)/momentum_ref),'-','linewidth',2, 'MarkerSize',10,'Color','#0072bd')

semilogy(t_j(1:nstep:end),abs((momentum_res_midpoint(1:nstep:end)-momentum_ref)/momentum_ref),'v','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
semilogy(t_j,abs((momentum_brunedschratz-momentum_ref)/momentum_ref),'-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e')
semilogy(t_j(1:nstep:end),abs((momentum_brunedschratz(1:nstep:end)-momentum_ref)/momentum_ref),'<','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
semilogy(t_j,abs((momentum_strang-momentum_ref)/momentum_ref),'-','linewidth',2, 'MarkerSize',10,'Color','#edb120')%'#8034eb')
semilogy(t_j(1:nstep:end),abs((momentum_strang(1:nstep:end)-momentum_ref)/momentum_ref),'o','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
semilogy(t_j,abs((momentum_lawson-momentum_ref)/momentum_ref),'-','linewidth',2, 'MarkerSize',10,'Color','black')%'#8034eb')
semilogy(t_j(1:nstep:end),abs((momentum_lawson(1:nstep:end)-momentum_ref)/momentum_ref),'d','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')


h(1)=semilogy(nan,nan,'v-','linewidth',2, 'MarkerSize',10,'Color','#0072bd','MarkerFaceColor','white')
h(2)=semilogy(nan,nan,'<-','linewidth',2, 'MarkerSize',10,'Color','#7e2f8e','MarkerFaceColor','white')
h(3)=semilogy(nan,nan,'o-','linewidth',2, 'MarkerSize',10,'Color','#edb120','MarkerFaceColor','white')
h(4)=semilogy(nan,nan,'d-','linewidth',2, 'MarkerSize',10,'Color','black','MarkerFaceColor','white')

hold off

tau
tau*N^2
set(gca,'FontSize',16)
%title('Strang, N=64')
xlabel('$t$','Interpreter','latex', 'FontSize', 22)
ylabel('Rel. error in momentum','Interpreter','latex', 'FontSize', 22)

legend(h,'Resonance-based midpoint','Bruned \& Schratz', 'Strang', 'Lawson','Interpreter','latex', 'FontSize', 16,'Position',[0.77 0.5 0.1 0.2])
grid on


ylim([1e-16,1e2])
hold off

set(gcf, 'Position',  [100, 100, 1400, 600])
%exportgraphics(gcf,strcat('images/kdv_momentum_preservation_',dataset,'.pdf'),'ContentType','vector')
