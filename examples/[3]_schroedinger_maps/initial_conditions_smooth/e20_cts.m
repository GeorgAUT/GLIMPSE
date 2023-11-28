function [out] = e20_cts(x)
%% Initial parallel frame e_2

out=sin(tau0_int_cts(x))*N0_cts(x)+cos(tau0_int_cts(x))*B0_cts(x);
end

