  %Generates state evolution and compares to observed MSEs for D-AMP and D-IT algorithms

addpath(genpath('..'));

%Parameters
denoiser='DnCNN';%DnCNN %fast-BM3D'
SamplingRate1=0.05 ;
SamplingRate2=0.1 ;
SamplingRate3=0.2 ;


%noise_sig=100;
SNR_dB=-5:5:30;
imsize=128;
filename='barbara.png';
iters=10;
N=300;%Number of tests to run.  Must be at least 2
n_DnCNN_layers=20;
LoadNetworkWeights(n_DnCNN_layers);% You need to run LoadNetworkWeights before you can use the DnCNN denoiser
% Im=double(imread(filename));
% x_0=imresize(Im,imsize/size(Im,1));

clean=test(64,64,1,3,[],[],12,12,20);% one user, three path, Az and ele random generate
[signal]=255*rescaleImage(clean);
%labels  = sigma_array/255.*randn(size(inputs),'single');%ÔëÉùµÄ´óÐ¡
% noisy=signal+sigma_true/255.*randn(size(clean),'single');
x_0=signal;

[height, width]=size(x_0);
n=length(x_0(:));
m1=round(n*SamplingRate1);
m2=round(n*SamplingRate2);
m3=round(n*SamplingRate3);


%Generate State Evolution and compute intermediate MSEs of D-AMP and D-IT
True_DAMP_MSE_array1=zeros(N,iters);
Predicted_MSE_array1=zeros(N,iters);

True_DAMP_MSE_array2=zeros(N,iters);
Predicted_MSE_array2=zeros(N,iters);

True_DAMP_MSE_array3=zeros(N,iters);
Predicted_MSE_array3=zeros(N,iters);

True_DAMP_MSE_array1(:,1)=mean(x_0(:).^2);
Predicted_MSE_array1(:,1)=mean(x_0(:).^2);

True_DAMP_MSE_array2(:,1)=mean(x_0(:).^2);
Predicted_MSE_array2(:,1)=mean(x_0(:).^2);

True_DAMP_MSE_array3(:,1)=mean(x_0(:).^2);
Predicted_MSE_array3(:,1)=mean(x_0(:).^2);



for jj=1:length(SNR_dB)
    jj
for ii=1:N
    
 SNR_linear=10.^(SNR_dB(jj)/10);

%     M=randn(m,n);
%     for j = 1:n
%         M(:,j) = M(:,j) ./ sqrt(sum(abs(M(:,j)).^2));
%     end
    
 M1 = (1/sqrt(m1))*(2*(rand(m1,n)<0.5) - 1);%
 M2 = (1/sqrt(m2))*(2*(rand(m2,n)<0.5) - 1);%
 M3 = (1/sqrt(m3))*(2*(rand(m3,n)<0.5) - 1);%
 
 
 % B=hadamard(n); %
%  M=B(randperm(n,m)',randperm(n));
    z1=M1*x_0(:);
    z2=M2*x_0(:);
    z3=M3*x_0(:);
    noise_power1=norm(z1)^2/length(z1)/SNR_linear;
    noise_power2=norm(z2)^2/length(z2)/SNR_linear;
    noise_power3=norm(z3)^2/length(z3)/SNR_linear;
    
    noise_sig1=sqrt(noise_power1)*randn(length(z1),1);
    noise_sig2=sqrt(noise_power2)*randn(length(z2),1);
    noise_sig3=sqrt(noise_power3)*randn(length(z3),1);
    
    
    
    y1=z1+noise_sig1;
    y2=z2+noise_sig2;
    y3=z3+noise_sig3;
    
    
    x1_t=zeros(n,1);
    
    x2_t=zeros(n,1);
    x3_t=zeros(n,1);
    
    z1_t=y1;
    z2_t=y2;
    z3_t=y3;
 
    for iter=2:iters
            [x1_tplus1,z1_tplus1,NA] = DAMP_oneIter(y1,z1_t,x1_t,width,height,denoiser,M1);
            z1_t=z1_tplus1; x1_t=x1_tplus1;
        True_DAMP_MSE_array1(ii,iter)=mean((x_0(:)-x1_tplus1(:)).^2);
     Predicted_MSE_array1(ii,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array1(ii,iter-1), m1,n,sqrt(noise_power1),denoiser,width,height);
     
     
      [x2_tplus1,z2_tplus1,NA] = DAMP_oneIter(y2,z2_t,x2_t,width,height,denoiser,M2);
            z2_t=z2_tplus1; x2_t=x2_tplus1;
        True_DAMP_MSE_array2(ii,iter)=mean((x_0(:)-x2_tplus1(:)).^2);
     Predicted_MSE_array2(ii,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array2(ii,iter-1), m2,n,sqrt(noise_power2),denoiser,width,height);
     
     
      [x3_tplus1,z3_tplus1,NA] = DAMP_oneIter(y3,z3_t,x3_t,width,height,denoiser,M3);
            z3_t=z3_tplus1; x3_t=x3_tplus1;
        True_DAMP_MSE_array3(ii,iter)=mean((x_0(:)-x3_tplus1(:)).^2);
     Predicted_MSE_array3(ii,iter)= DAMP_SE_Prediction(x_0(:), Predicted_MSE_array3(ii,iter-1), m3,n,sqrt(noise_power3),denoiser,width,height);
     
     
     
end
end
MSE_low_True1=min(True_DAMP_MSE_array1,[],2);
MSE_low_Predict1=min(Predicted_MSE_array1,[],2);
True_DAMP_MSE_snr1(jj)=sum(MSE_low_True1)/N/mean((x_0(:)))^2;
Predicted_MSE_snr1(jj)=sum(MSE_low_Predict1)/N/mean((x_0(:)))^2;


MSE_low_True2=min(True_DAMP_MSE_array2,[],2);
MSE_low_Predict2=min(Predicted_MSE_array2,[],2);
True_DAMP_MSE_snr2(jj)=sum(MSE_low_True2)/N/mean((x_0(:)))^2;
Predicted_MSE_snr2(jj)=sum(MSE_low_Predict2)/N/mean((x_0(:)))^2;

MSE_low_True3=min(True_DAMP_MSE_array3,[],2);
MSE_low_Predict3=min(Predicted_MSE_array3,[],2);
True_DAMP_MSE_snr3(jj)=sum(MSE_low_True3)/N/mean((x_0(:)))^2;
Predicted_MSE_snr3(jj)=sum(MSE_low_Predict3)/N/mean((x_0(:)))^2;




end
True_DAMP_MSE1=mean(True_DAMP_MSE_array1)';
Predicted_MSE1=mean(Predicted_MSE_array1)';

True_DAMP_MSE2=mean(True_DAMP_MSE_array2)';
Predicted_MSE2=mean(Predicted_MSE_array2)';

True_DAMP_MSE3=mean(True_DAMP_MSE_array3)';
Predicted_MSE3=mean(Predicted_MSE_array3)';




% True_DAMP_MSE_std=std(True_DAMP_MSE_array)';
% True_DIT_MSE=mean(True_DIT_MSE_array)';
% True_DIT_MSE_std=std(True_DIT_MSE_array)';

% Predicted_MSE_std=std(Predicted_MSE_array)';




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
plot(SNR_dB,10*log10(True_DAMP_MSE_snr1),'ro');hold on 
plot(SNR_dB,10*log10(Predicted_MSE_snr1),'r-');hold on

plot(SNR_dB,10*log10(True_DAMP_MSE_snr2),'bo');hold on 
plot(SNR_dB,10*log10(Predicted_MSE_snr2),'b-');hold on

plot(SNR_dB,10*log10(True_DAMP_MSE_snr3),'ko');hold on 
plot(SNR_dB,10*log10(Predicted_MSE_snr3),'k-');hold on
%figure(2)
% plot(0:29,10*log10(True_DAMP_MSE),'ro');hold on 
% plot(0:29,10*log10(Predicted_MSE),'r-');hold on

grid on
xlabel('SNR (dB)');
ylabel('NMSE (dB)');hold on