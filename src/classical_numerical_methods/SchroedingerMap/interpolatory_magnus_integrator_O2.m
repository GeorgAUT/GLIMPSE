function [out] = interpolatory_magnus_integrator_O2(y,u_hat_old,u_hat,h)
%Interpolatory Magnus integrator of order 2 as introduced for SM equation
%in arXiv:2211.01282
%   y=[T0;e10;e20];

N=max(size(y));

M0=zeros(3,3,N);
M1=zeros(3,3,N);

dxu_old=ifftpi(dx(u_hat_old));
dxu=ifftpi(dx(u_hat));

M0(1,2,:)=-imag(dxu_old);
M0(1,3,:)=real(dxu_old);
M0(2,1,:)=imag(dxu_old);
M0(3,1,:)=-real(dxu_old);
M0(2,3,:)=-abs(ifftpi(u_hat_old)).^2/2;
M0(3,2,:)=abs(ifftpi(u_hat_old)).^2/2;

M1(1,2,:)=-imag(dxu);
M1(1,3,:)=real(dxu);
M1(2,1,:)=imag(dxu);
M1(3,1,:)=-real(dxu);
M1(2,3,:)=-abs(ifftpi(u_hat)).^2/2;
M1(3,2,:)=abs(ifftpi(u_hat)).^2/2;

Mat_t1=h/2*(M0+M1);
Mat_t2=h^2/6*(pagemtimes(M1,M0)-pagemtimes(M0,M1));

out=pagemtimes(pagemexpm(Mat_t1+Mat_t2),y);

end

