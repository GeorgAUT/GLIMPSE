clear all;

load('data/numerical_experiments_H1_data.mat') % Load data from experiments
%% Simple convergence rate plots
figure(1)

loglog(tau_j, error_first_order)
hold on
loglog(tau_j, error_brunedschratz)
loglog(tau_j, error_symmetric)
loglog(tau_j, error_lie)
loglog(tau_j, error_strang)
ylabel('L^2-error')
xlabel('tau')
legend('Ostermann & Schratz', 'Bruned & Schratz', 'Symmetrised Res.-Based', 'Lie', 'Strang')
grid on
hold off

exportgraphics(gcf,'images/convergence_rates_example.pdf','ContentType','vector')

%% Simple CPU-time plots
figure(2)
semilogy(cputime_first_order, error_first_order)
hold on
semilogy(cputime_brunedschratz, error_brunedschratz)
semilogy(cputime_symmetric, error_symmetric)
semilogy(cputime_lie, error_lie)
semilogy(cputime_strang, error_strang)
ylabel('L^2-error')
xlabel('CPU-time (sec)')
legend('Ostermann & Schratz', 'Bruned & Schratz', 'Symmetrised Res.-Based', 'Lie', 'Strang')
grid on
hold off

exportgraphics(gcf,'images/cput_times_example.pdf','ContentType','vector')
