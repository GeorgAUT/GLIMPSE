function [out] = tau0_cts(x)
%% Initial condition for tau (torsion)
out=2*cos(x).*(1 + 1./(3 - 2*cos(2*x)));
end

