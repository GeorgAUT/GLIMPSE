function [out] = u0_cts(x)
%% Initial condition for u0 (NLS flow)
out= kappa0_cts(x).*exp(i*tau0_int_cts(x));
end

