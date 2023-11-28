function [out] = SME_low_reg_integrator(y,u_hat_t1,u_hat_t2,h)
%% Take a timestep with the FLowRH transform as introduced in arXiv:2211.01282
% Note this code does not perform the fast version of the algorithm
%
%   Input:  y=[T0;e10;e20]....value of y at time tn
%           u_hat_t1....value of NLS solution at time tn
%           u_hat_t2....value of NLS solution at time tn+h
%           h....time step

N=max(size(y));

% Highly oscillatory integrator for the easy terms:
dxv_mixed_hat=dx(u_hat_t1+expilaplacian(u_hat_t2,-h))/2;
dxu_int=ifftpi(h*phi1laplacian(dxv_mixed_hat,h));

% Highly oscillatory integrator for the hard terms:
v_mixed_hat=(expilaplacian(u_hat_t1,-h)+u_hat_t2)/2;
u_mixed_hat=(u_hat_t1+u_hat_t2)/2;

v_mixed_conj_hat=fftpi(conj(ifftpi(v_mixed_hat)));
u_mixed_conj_hat=fftpi(conj(ifftpi(u_mixed_hat)));

int_out=zeros(size(u_hat_t1));
int_out2=zeros(size(u_hat_t1));


parfor l=1:N

    k1=transpose(max(1,l-N/2):min(N,l+N/2-1));

    k2=l-k1;
    k1=k1-N/2;

    a1=(exp(-i*h*(k1.^2))-1)./(-i*(k1.^2)).*v_mixed_hat(k1+N/2).*u_mixed_conj_hat(k2+N/2);
    a10=(exp(-i*h*(k1.^2-k2.^2))-1)./(-i*(k1.^2-k2.^2)).*v_mixed_hat(k1+N/2).*v_mixed_conj_hat(k2+N/2);

    a2=(exp(-i*h*(-k2.^2))-1)./(-i*(-k2.^2)).*u_mixed_hat(k1+N/2).*v_mixed_conj_hat(k2+N/2);
    a20=(exp(-i*h*(k1.^2-k2.^2))-1)./(-i*(k1.^2-k2.^2)).*v_mixed_hat(k1+N/2).*v_mixed_conj_hat(k2+N/2);
    
    b=h*v_mixed_hat(k1+N/2).*v_mixed_conj_hat(k2+N/2);

    a=zeros(size(a1));
    a(find((k1.^2)>=k2.^2))=a1(find((k1.^2)>=k2.^2));
    a(find((k1.^2)==0))=b(find((k1.^2)==0));
    a(find((k1.^2)<k2.^2))=a2(find((k1.^2)<k2.^2));
    int_out(l)=sum(a);

end

int_out=ifftpi(int_out)/2;


Mat_t1=zeros(3,3,N);
Mat_t1(1,2,:)=-imag(dxu_int);
Mat_t1(1,3,:)=real(dxu_int);
Mat_t1(2,1,:)=imag(dxu_int);
Mat_t1(3,1,:)=-real(dxu_int);
Mat_t1(2,3,:)=-int_out;
Mat_t1(3,2,:)=int_out;

out=pagemtimes(pagemexpm(Mat_t1),y);

end

