function [output]=pagecross(A,B)
% Pagecross computes coordinatewise cross product of row-vectors
N=size(A,3);
if size(A)~=size(B)
    fprintf('Missmatching sizes in pagecross')
    return
end
output=zeros(size(A));
parfor j=1:N
    output(:,:,j)=cross(A(:,:,j),B(:,:,j));
end
end