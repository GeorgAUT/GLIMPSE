function [out]=fourth_order_splitting_NLS(u,h)
%% Compute fourth order splitting
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}

a1=1.3512071919596577718181;
a2=1-2*a1;

out=strang_NLS(strang_NLS(strang_NLS(u,a1*h),a2*h),a1*h);
        
end