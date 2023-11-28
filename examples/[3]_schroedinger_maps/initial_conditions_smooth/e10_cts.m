function [out] = e10_cts(x)
%% Initial parallel frame e_1

out= cos(tau0_int_cts(x))*N0_cts(x)-sin(tau0_int_cts(x))*B0_cts(x);
end

