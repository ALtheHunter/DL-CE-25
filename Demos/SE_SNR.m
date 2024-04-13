  %Generates state evolution and compares to observed MSEs for D-AMP and D-IT algorithms

addpath(genpath('..'));

%Parameters
denoiser='DnCNN';%DnCNN %fast-BM3D'
SamplingRate=0.05 ;
%noise_sig=100;
SNR_dB=0:5:30;
imsize=128;
filename='barbara.png';
iters=10;
N=1000;%Number of tests to run.  Must be at least 2
n_DnCNN_layers=20;
LoadNetworkWeights(n_DnCNN_layers);% You need to run LoadNetworkWeights before you can use the DnCNN denoiser

clean=test(64,64,1,4,[],[],12,12,20);% one user, three path, Az and ele random generate
[signal]=255*rescaleImage(clean);

x_0=signal;

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
for jj=1:length(SNR_dB)
    jj
for ii=1:N
    
 SNR_linear=10.^(SNR_dB(jj)/10);
 M = (1/sqrt(m))*(2*(rand(m,n)<0.5) - 1);%
     z=M*x_0(:);
    noise_power=norm(z)^2/length(z)/SNR_linear;
    noise_sig=sqrt(noise_power)*randn(length(z),1);
    y=z+noise_sig;
    x_t=zeros(n,1);x_t2=zeros(n,1);
    z_t=y;
    z_t2=y;
    for iter=2:iters
            [x_tplus1,z_tplus1,NA] = DAMP_oneIter(y,z_t,x_t,width,height,denoiser,M);
            z_t=z_tplus1; x_t=x_tplus1;
            True_DAMP_MSE_array(ii,iter)=mean((x_0(:)-x_tplus1(:)).^2);
            Predicted_MSE_array(ii,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array(ii,iter-1), m,n,sqrt(noise_power),denoiser,width,height);
    end
end
MSE_low_True=min(True_DAMP_MSE_array,[],2);
MSE_low_Predict=min(Predicted_MSE_array,[],2);
True_DAMP_MSE_snr(jj)=sum(MSE_low_True)/N/mean((x_0(:)))^2;

Predicted_MSE_snr(jj)=sum(MSE_low_Predict)/N/mean((x_0(:)))^2;
end

plot(SNR_dB,10*log10(True_DAMP_MSE_snr),'ro');hold on 
plot(SNR_dB,10*log10(Predicted_MSE_snr),'r-');hold on

%%
%Parameters
denoiser='DnCNN';%DnCNN %fast-BM3D'
SamplingRate=0.1 ;
n_DnCNN_layers=20;
LoadNetworkWeights(n_DnCNN_layers);% You need to run LoadNetworkWeights before you can use the DnCNN denoiser

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
for jj=1:length(SNR_dB)
    jj
for ii=1:N
    
 SNR_linear=10.^(SNR_dB(jj)/10);
 M = (1/sqrt(m))*(2*(rand(m,n)<0.5) - 1);%
     z=M*x_0(:);
    noise_power=norm(z)^2/length(z)/SNR_linear;
    noise_sig=sqrt(noise_power)*randn(length(z),1);
    y=z+noise_sig;
    x_t=zeros(n,1);x_t2=zeros(n,1);
    z_t=y;
    z_t2=y;
    for iter=2:iters
            [x_tplus1,z_tplus1,NA] = DAMP_oneIter(y,z_t,x_t,width,height,denoiser,M);
            z_t=z_tplus1; x_t=x_tplus1;
            True_DAMP_MSE_array(ii,iter)=mean((x_0(:)-x_tplus1(:)).^2);
            Predicted_MSE_array(ii,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array(ii,iter-1), m,n,sqrt(noise_power),denoiser,width,height);
    end
end
MSE_low_True=min(True_DAMP_MSE_array,[],2);
MSE_low_Predict=min(Predicted_MSE_array,[],2);
True_DAMP_MSE_snr(jj)=sum(MSE_low_True)/N/mean((x_0(:)))^2;

Predicted_MSE_snr(jj)=sum(MSE_low_Predict)/N/mean((x_0(:)))^2;
end

plot(SNR_dB,10*log10(True_DAMP_MSE_snr),'bo');hold on 
plot(SNR_dB,10*log10(Predicted_MSE_snr),'b-');hold on

%%%%
%Parameters
denoiser='DnCNN';%DnCNN %fast-BM3D'
SamplingRate=0.2 ;
n_DnCNN_layers=20;
LoadNetworkWeights(n_DnCNN_layers);% You need to run LoadNetworkWeights before you can use the DnCNN denoiser

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
for jj=1:length(SNR_dB)
    jj
for ii=1:N
    
 SNR_linear=10.^(SNR_dB(jj)/10);
 M = (1/sqrt(m))*(2*(rand(m,n)<0.5) - 1);%
     z=M*x_0(:);
    noise_power=norm(z)^2/length(z)/SNR_linear;
    noise_sig=sqrt(noise_power)*randn(length(z),1);
    y=z+noise_sig;
    x_t=zeros(n,1);x_t2=zeros(n,1);
    z_t=y;
    z_t2=y;
    for iter=2:iters
            [x_tplus1,z_tplus1,NA] = DAMP_oneIter(y,z_t,x_t,width,height,denoiser,M);
            z_t=z_tplus1; x_t=x_tplus1;
            True_DAMP_MSE_array(ii,iter)=mean((x_0(:)-x_tplus1(:)).^2);
            Predicted_MSE_array(ii,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array(ii,iter-1), m,n,sqrt(noise_power),denoiser,width,height);
    end
end
MSE_low_True=min(True_DAMP_MSE_array,[],2);
MSE_low_Predict=min(Predicted_MSE_array,[],2);
True_DAMP_MSE_snr(jj)=sum(MSE_low_True)/N/mean((x_0(:)))^2;

Predicted_MSE_snr(jj)=sum(MSE_low_Predict)/N/mean((x_0(:)))^2;
end

plot(SNR_dB,10*log10(True_DAMP_MSE_snr),'ko');hold on 
plot(SNR_dB,10*log10(Predicted_MSE_snr),'k-');hold on
