%Demonstrates compressively sampling and LD(V)AMP recovery of an image.
%Requires  matconvnet and gampmatlab in the path
addpath(genpath('..'));
% addpath(genpath('~/gampmatlab'));
% addpath('~/matconvnet/matlab');
clear all
%Parameters
denoiser1='BM3D';%Available options are NLM, Gauss, Bilateral, BLS-GSM, BM3D, fast-BM3D, BM3D-SAPCA, and DnCNN
denoiser2='DnCNN';
filename='barbara.png';
SamplingRate=0.05;
AMP_iters=20;
VAMP_iters=20;
imsize=128; %512
n_DnCNN_layers=17;%Other option is 17
LoadNetworkWeights(n_DnCNN_layers);
SNR_dB=5:1:20;


num=100;
for jj=1:length(SNR_dB)
    jj
sigma_noise_de=10^(SNR_dB(jj)/10);
for ii=1:num
    
h_lens=test(64,64,1,3,[],[],12,12,20);% one user, three path, Az and ele random generate
%[x_0]=h_lens;
[x_0]=255*rescaleImage(h_lens);

% ImIn=double(imread(filename));
% x_0=imresize(ImIn,imsize/size(ImIn,1));
% [x_0]=rescaleImage(x_0);

[height, width]=size(x_0);
n=length(x_0(:));
m=round(n*SamplingRate);
errfxn = @(x_hat) PSNR(x_0,reshape(x_hat,[height width]));
MSEfxn = @(x_hat) MSE_cal(x_0,reshape(x_hat,[height width]));

%Generate Gaussian Measurement Matrix
  A=1/sqrt(m)*randn(m,n);

%  B=hadamard(n); %
%  A=1/sqrt(m)*B(randperm(n,m)',randperm(n));

M=@(x) A*x;
Mt=@(x) A.'*x;

%Generate Coded Diffraction Pattern Measurement Matrix
% signvec = exp(1i*2*pi*rand(n,1));
% inds=[1;randsample(n-1,m-1)+1];
% I=speye(n);
% SubsampM=I(inds,:);
% M=@(x) SubsampM*reshape(fft2(reshape(bsxfun(@times,signvec,x(:)),[height,width])),[n,1])*(1/sqrt(n))*sqrt(n/m);
% Mt=@(x) bsxfun(@times,conj(signvec),reshape(ifft2(reshape(SubsampM'*x(:),[height,width])),[n,1]))*sqrt(n)*sqrt(n/m);
U=@(x) x(:);
Ut= @(x) x(:);
d=ones(m,1)*n/m;

%Compressively sample the image
z=M(x_0(:));% noiseless observation
%noise=(norm(z)/length(z)/sigma_noise_de/2)*(randn(length(z),1)+j*randn(length(z),1));
noise=sqrt((norm(z)^2/length(z)/sigma_noise_de))*randn(length(z),1);
y=z+noise;
%y=z;
%Recover Signal using D-AMP algorithms
t0=tic;[x_hat1,mmse1(jj,ii)]  = DAMP_SNR(y,AMP_iters,height,width,denoiser1,M,Mt,errfxn,MSEfxn);t1=toc(t0);
% t0=tic;[x_hat2,psnr2(ii,:)] = DVAMP(y,VAMP_iters,height,width,denoiser1,M,Mt,errfxn, U, Ut, d);t2=toc(t0);

%D(V)AMP Recovery Performance
performance1=PSNR(x_0,x_hat1);
% performance2=PSNR(x_0,x_hat2);

%MSE_DAMP(ii,:)
end
end
% display([num2str(SamplingRate*100),'% Sampling ', denoiser1, '-AMP: PSNR=',num2str(performance1),', time=',num2str(t1)])
% display([num2str(SamplingRate*100),'% Sampling ', denoiser2, '-VAMP: PSNR=',num2str(performance2),', time=',num2str(t2)])

% psnr1_ave=sum(psnr1,1)/num;
mmse1_ave=sum(mmse1,2)/num;
% psnr2_ave=sum(psnr2,1)/num;
%Plot Recovered Signals
% figure(1); clf;
% subplot(1,3,1);
% imshow(uint8(x_0));title('Original Image');
% subplot(1,3,2);
% imshow(uint8(x_hat1));title([denoiser1, '-AMP']);
% subplot(1,3,3);
% imshow(uint8(x_hat2));title([denoiser2, '-VAMP']);
% 
% figure(1); clf;
% subplot(1,3,1);
% imshow(x_0,[]);title('Original Image');
% subplot(1,3,2);
% imshow(x_hat1,[]);title([denoiser1, '-AMP']);
% subplot(1,3,3);
% imshow(x_hat2,[]);title([denoiser2, '-VAMP']);

%Plot PSNR Trajectories
%figure(1); clf;
% plot(psnr1_ave,'.-','Displayname',[denoiser1,'-AMP (20)'])
% hold on; 

plot(SNR_dB,10*log10(mmse1_ave),'.-','Displayname',[denoiser1,'-AMP-',num2str(n_DnCNN_layers)])
hold on; 
% plot(psnr2_ave,'.-','Displayname',[denoiser2,'-VAMP (20)']); 
% hold on;
grid on;
legend(gca,'Show','Location','SouthEast')
xlabel('iteration')
ylabel('NMSE')