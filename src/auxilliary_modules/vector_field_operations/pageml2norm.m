function [norm] = pageml2norm(T)
% Computes L2 norm of vector field
norm=(pagemtimes(T,pagetranspose(T))).^(1/2);
end