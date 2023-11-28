function [out] = interpolatory_magnus_integrator_O4(y,u_hat_t1,u_hat_t2,h)
%Interpolatory Magnus integrator of order 4 as introduced for SM equation
%in arXiv:2211.01282
%   y=[T0;e10;e20];

N=max(size(y));

M1=zeros(3,3,N);
M2=zeros(3,3,N);

dxu_t1=ifftpi(dx(u_hat_t1));
dxu_t2=ifftpi(dx(u_hat_t2));

M1(1,2,:)=-imag(dxu_t1);
M1(1,3,:)=real(dxu_t1);
M1(2,1,:)=imag(dxu_t1);
M1(3,1,:)=-real(dxu_t1);
M1(2,3,:)=-abs(ifftpi(u_hat_t1)).^2/2;
M1(3,2,:)=abs(ifftpi(u_hat_t1)).^2/2;

M2(1,2,:)=-imag(dxu_t2);
M2(1,3,:)=real(dxu_t2);
M2(2,1,:)=imag(dxu_t2);
M2(3,1,:)=-real(dxu_t2);
M2(2,3,:)=-abs(ifftpi(u_hat_t2)).^2/2;
M2(3,2,:)=abs(ifftpi(u_hat_t2)).^2/2;

% weights
w1=0.5;
w2=0.5;
w11=0.125;
w12=-0.0193375672974064397036;
w21=0.2693375672974066548093;
w22=0.125;
w111=0.0208333333333333599324;
w112=-0.0027510584910182774888;
w121=-0.0041666666666666666088;
w122=-0.0027510584910182735857;
w211=0.0694177251576849480008;
w212=-0.0041666666666666640068;
w221=0.0694177251576849896342;
w222=0.0208333333333333495241;

Mat_t1=h*(w1*M1+w2*M2);
Mat_t2=h^2/2*(w12-w21)*pagecommutator(M1,M2);
% Mat_t3=h^3/6*w112*(pagecommutator(M1,pagecommutator(M1,M2)))...
%     +h^3/6*w121*(pagecommutator(M1,pagecommutator(M2,M1))-pagecommutator(M1,pagecommutator(M2,M1)))...
%     +h^3/6*w122*(-pagecommutator(M2,pagecommutator(M2,M1)))...
%     +h^3/6*w211*(-pagecommutator(M1,pagecommutator(M1,M2)))...
%     +h^3/6*w212*(pagecommutator(M2,pagecommutator(M1,M2))-pagecommutator(M2,pagecommutator(M1,M2)))...
%     +h^3/6*w221*(pagecommutator(M2,pagecommutator(M2,M1)));
% Mat_t3=h^3/6*w111*(pagecommutator(M1,pagecommutator(M1,M1))-pagecommutator(M1,pagecommutator(M1,M1)))...
%     +h^3/6*w112*(pagecommutator(M1,pagecommutator(M1,M2))-pagecommutator(M2,pagecommutator(M1,M1)))...
%     +h^3/6*w121*(pagecommutator(M1,pagecommutator(M2,M1))-pagecommutator(M1,pagecommutator(M2,M1)))...
%     +h^3/6*w122*(pagecommutator(M1,pagecommutator(M2,M2))-pagecommutator(M2,pagecommutator(M2,M1)))...
%     +h^3/6*w211*(pagecommutator(M2,pagecommutator(M1,M1))-pagecommutator(M1,pagecommutator(M1,M2)))...
%     +h^3/6*w212*(pagecommutator(M2,pagecommutator(M1,M2))-pagecommutator(M2,pagecommutator(M1,M2)))...
%     +h^3/6*w221*(pagecommutator(M2,pagecommutator(M2,M1))-pagecommutator(M1,pagecommutator(M2,M2)))...
%     +h^3/6*w222*(pagecommutator(M2,pagecommutator(M2,M2))-pagecommutator(M2,pagecommutator(M2,M2)));

out=pagemtimes(pagemexpm(Mat_t1+Mat_t2),y);

end

