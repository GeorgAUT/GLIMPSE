function [out] = N0_cts(x)
%% Initial normal vector field
out=[-cos(x)+3*cos(3*x),-sin(x)+3*sin(3*x),-2*sin(x)]./2./sqrt(3-2*cos(2*x));
end

