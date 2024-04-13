%%% Test the model
addpath('~/matconvnet/matlab');
addpath('../../../Utils')
addpath('../../../TestImages')
vl_setupnn

%Parameters
sigma_true=120;
sigma_min = 100;
sigma_max = 150;
showResult  = 1;
useGPU      = 1;
pauseTime   = 0;

%Select which network to use
modelName = ['DnCNN_26Layer_sigma',num2str(sigma_min),'to',num2str(sigma_max)]; %%% model name
epoch = 40;

%%% load Gaussian denoising model
load(fullfile('NewNetworks',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
%load(fullfile('NewNetworks',[modelName,'-best','.mat']));
net = vl_simplenn_tidy(net);
net.layers = net.layers(1:end-1);

%%%
net = vl_simplenn_tidy(net);

%%% move to gpu
if useGPU
    net = vl_simplenn_move(net, 'gpu') ;
end

%Create the noisy image
%clean=double(imread('barbara.png'));

% noisy=clean+40*randn(size(clean));
% noisy=single(noisy);

clean=test(50,50,1,3,[],[],12,12,20);% one user, three path, Az and ele random generate
[signal]=255*rescaleImage(clean);
%labels  = sigma_array/255.*randn(size(inputs),'single');%ÔëÉùµÄ´óÐ¡
noisy=signal+sigma_true/255.*randn(size(clean),'single');
noisy=single(noisy);

%%% convert to GPU
if useGPU
    noisy = gpuArray(noisy);
end

res    = vl_simplenn(net,noisy,[],[],'conserveMemory',true,'mode','test');
denoised = noisy - res(end).x;

%%% convert to CPU
if useGPU
    denoised = gather(denoised);
    noisy  = gather(noisy);
end

%%% calculate PSNR and SSIM
Input_PSNR=PSNR(denoised,signal);
Output_PSNR=PSNR(noisy,signal);

figure(1);title('DnCNN-new')
subplot(1,3,1);imshow(signal,[]);title('signal');
subplot(1,3,2);imshow(noisy,[]);title('Noisy input');
subplot(1,3,3);imshow(denoised,[]);title('Denoised output');

% display('Input psnr: ',Input_PSNR);
%  display('Output psnr: ',Output_PSNR);

disp(Input_PSNR);
disp(Output_PSNR);

%%
% modelName = ['DnCNN_20Layer_sigma',num2str(sigma_min),'to',num2str(sigma_max)]; %%% model name
% epoch = 50;
% %%% load Gaussian denoising model
% %load(fullfile('NewNetworks',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
% %load(fullfile('NewNetworks',[modelName,'-best','.mat']));
% load(fullfile('NewNetworks',modelName,[modelName,'-epoch-',num2str(epoch),'.mat']));
% net = vl_simplenn_tidy(net);
% net.layers = net.layers(1:end-1);
% 
% %%%
% net = vl_simplenn_tidy(net);
% 
% %%% move to gpu
% if useGPU
%     net = vl_simplenn_move(net, 'gpu') ;
% end
% 
% 
% 
% %%% convert to GPU
% if useGPU
%     noisy = gpuArray(noisy);
% end
% 
% res    = vl_simplenn(net,noisy,[],[],'conserveMemory',true,'mode','test');
% denoised = noisy - res(end).x;
% 
% %%% convert to CPU
% if useGPU
%     denoised = gather(denoised);
%     noisy  = gather(noisy);
% end
% 
% %%% calculate PSNR and SSIM
% Input_PSNR=PSNR(denoised,clean);
% Output_PSNR=PSNR(noisy,clean);
% 
% figure(2);title('DnCNN-old')
% subplot(1,3,1);imshow(signal,[]);title('signal');
% subplot(1,3,2);imshow(noisy,[]);title('Noisy input');
% subplot(1,3,3);imshow(denoised,[]);title('Denoised output');
% 
% % display('Input psnr: ',Input_PSNR);
% %  display('Output psnr: ',Output_PSNR);
% 
% disp(Input_PSNR);
% disp(Output_PSNR);
