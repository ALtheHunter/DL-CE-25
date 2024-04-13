%Generates state evolution and compares to observed MSEs for D-AMP and D-IT algorithms

addpath(genpath('..'));

%Parameters
denoiser='DnCNN';%DnCNN %fast-BM3D'
SamplingRate=0.05;
SNR_dB=10;
SNR=10^(SNR_dB/10);
imsize=128;
filename='barbara.png';
iters=10;
N=400;%Number of tests to run.  Must be at least 2
n_DnCNN_layers=20;
LoadNetworkWeights(n_DnCNN_layers);% You need to run LoadNetworkWeights before you can use the DnCNN denoiser
% Im=double(imread(filename));
% x_0=imresize(Im,imsize/size(Im,1));

clean=test(64,64,1,4,[],[],12,12,20);% one user, three path, Az and ele random generate
[signal]=rescaleImage(clean);
%labels  = sigma_array/255.*randn(size(inputs),'single');%ÔëÉùµÄ´óÐ¡
% noisy=signal+sigma_true/255.*randn(size(clean),'single');
x_0=255*signal;

[height, width]=size(x_0);
n=length(x_0(:));
m=round(n*SamplingRate);

%Generate State Evolution and compute intermediate MSEs of D-AMP and D-IT
True_DAMP_MSE_array=zeros(N,iters);
True_DIT_MSE_array=zeros(N,iters);
Predicted_MSE_array=zeros(N,iters);
True_DAMP_MSE_array(:,1)=mean(x_0(:).^2);
True_DIT_MSE_array(:,1)=mean(x_0(:).^2);
Predicted_MSE_array(:,1)=mean(x_0(:).^2);
for i=1:N
    i
%     M=randn(m,n);
%     for j = 1:n
%         M(:,j) = M(:,j) ./ sqrt(sum(abs(M(:,j)).^2));
%     end
% B=hadamard(n); %
%  M=B(randperm(n,m)',randperm(n));

    M = (1/sqrt(m))*(2*(rand(m,n)<0.5) - 1);%
 % B=hadamard(n); %
%  M=B(randperm(n,m)',randperm(n));
    z=M*x_0(:);
    noise_power=norm(z)^2/length(z)/SNR;
    noise_sig=sqrt(noise_power);



    y=z+noise_sig*randn(m,1);
    x_t=zeros(n,1);x_t2=zeros(n,1);
    z_t=y;
    z_t2=y;
    for iter=2:iters
            [x_tplus1,z_tplus1,NA] = DAMP_oneIter(y,z_t,x_t,width,height,denoiser,M);
            z_t=z_tplus1; x_t=x_tplus1;
            [x_tplus12,z_tplus12,NA] = DIT_oneIter(y,z_t2,x_t2,width,height,denoiser,M);
            z_t2=z_tplus12; x_t2=x_tplus12;
            True_DAMP_MSE_array(i,iter)=mean((x_0(:)-x_tplus1(:)).^2);
            True_DIT_MSE_array(i,iter)=mean((x_0(:)-x_tplus12(:)).^2);
            Predicted_MSE_array(i,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array(i,iter-1), m,n,noise_sig,denoiser,width,height);
    end
end
True_DAMP_MSE=mean(True_DAMP_MSE_array)'/mean(x_0(:).^2);
True_DAMP_MSE_std=std(True_DAMP_MSE_array)';
True_DIT_MSE=mean(True_DIT_MSE_array)';
True_DIT_MSE_std=std(True_DIT_MSE_array)';
Predicted_MSE=mean(Predicted_MSE_array)'/mean(x_0(:).^2);
Predicted_MSE_std=std(Predicted_MSE_array)';

%Plot Results
% h=figure; 
% hold
% errorbar(0:29,Predicted_MSE,Predicted_MSE_std,'-.b');
% errorbar(0:29,True_DAMP_MSE,True_DAMP_MSE_std,'-g');
% errorbar(0:29,True_DIT_MSE,True_DIT_MSE_std,'-r');
% title([denoiser,'-AMP and ', denoiser '-IT State Evolution']);
% legend(['Predicted ',num2str(100*SamplingRate),'%'], [denoiser,'-AMP ', num2str(100*SamplingRate),'%'], [denoiser,'-IT ', num2str(100*SamplingRate),'%']);
% xlabel('Iteration');
% ylabel('MSE');

%figure(2)
% semilogy(0:29,True_DAMP_MSE,'bo');hold on 
% semilogy(0:29,Predicted_MSE,'bx-');hold on

plot(1:10,10*log10(True_DAMP_MSE),'ro');hold on 
plot(1:10,10*log10(Predicted_MSE),'rx-');hold on
grid on
xlabel('iteration');
ylabel('NMSE (dB)');hold on