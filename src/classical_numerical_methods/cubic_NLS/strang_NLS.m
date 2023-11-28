function [out]=strang_NLS(u,h)
%% Compute strang splitting
% Input:    u...(\hat{u}_n)_{n=-N/2+1}^{N/2}

out=expilaplacian(nls_nonlinear_part_SME(expilaplacian(u,h/2),h),h/2);
        
end