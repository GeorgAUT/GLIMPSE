function [out] = T0_cts(x)
%% Initial tangent vector field
out=[cos(2*x).*sin(x),sin(2*x).*sin(x),cos(x)];
end

