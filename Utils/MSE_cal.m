function [MSE]=MSE_cal(x,x_hat)
% function [psnr]=PSNR(x,x_hat)
% PSNR computes the psnr of x_hat; an estimate of x.
% Input:
%       x     : the true signal
%       x_hat : an estimate of x
%Output:
%       psnr  : the PSNR of x_hat
    [imheight, imwidth]=size(x);
  MSE=norm(x(:)-x_hat(:))^2/norm(x(:))^2;
  