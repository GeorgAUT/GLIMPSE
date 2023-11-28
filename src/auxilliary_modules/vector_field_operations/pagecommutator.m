function [out] = pagecommutator(M1,M0)
% Computes commutator of two matrixfields/frames
out=pagemtimes(M1,M0)-pagemtimes(M0,M1);
end

