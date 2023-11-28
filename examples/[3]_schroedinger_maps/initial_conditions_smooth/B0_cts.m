function [out] = B0_cts(x)
%% Initial condition for binormal vector field

out=[-8*cos(x).^3.*sin(x),-1+2*cos(2*x)+cos(4*x),4*sin(x).^2]./2./sqrt(3-2*cos(2*x));
end

