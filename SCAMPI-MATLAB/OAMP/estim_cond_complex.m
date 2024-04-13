function [umean,uvar]=estim_cond_complex(Xmod,cmean,cvar)
%cmean:x_AB,vector
%vmean:v_AB,scalar
%Xmod:xn
     cvar = max(cvar, eps) ;     
     [size_n,size_m] = size(cmean) ; 
     cmean = reshape(cmean, size_m*size_n, 1) ;
     %---------------------------------------------------
     logpxr = bsxfun(@times, -1./cvar, abs(bsxfun(@minus, cmean, Xmod)).^2);
     pxr = exp(logpxr);
     z=sum(pxr,2);
     z(z==0)=eps;
     pxr = bsxfun(@rdivide, pxr,z);
     umean = sum(bsxfun(@times, pxr, Xmod),2);  
     uvar = sum(pxr .* abs(bsxfun(@minus, umean, Xmod)) .^2, 2);
     uvar = sum(uvar)/length(uvar);
end