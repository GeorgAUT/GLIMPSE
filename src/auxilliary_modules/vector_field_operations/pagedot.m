function [output]=pagedot(A,B)
% Pagecross computes coordinatewise dot product of row-vectors
N=size(A,3);
if size(A)~=size(B)
    fprintf('Missmatching sizes in pagedot')
    return
end
output=zeros(size(A));
parfor j=1:N
    output(:,:,j)=dot(A(:,:,j),B(:,:,j));
end
end