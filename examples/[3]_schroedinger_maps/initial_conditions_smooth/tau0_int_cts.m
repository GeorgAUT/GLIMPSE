function [out] = tau0_int_cts(x)
%% Initial condition for \int tau
out= (atan(2*sin(x)) + 2*sin(x));
end

