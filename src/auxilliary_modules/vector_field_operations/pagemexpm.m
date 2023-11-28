function output=pagemexpm(Mat)
%% Matrix exponential in parallel along third coordinate
L=size(Mat,3);
output=zeros(size(Mat));
parfor l=1:L
    output(:,:,l)=expm(Mat(:,:,l));
end
end
