clear all 
clc
had = 1; % Do you want to use an Hadamard or a random measurement matrix?
f = 2;
N = (32*f)^2; % size of image (power of two with Hadamard operator)
% image_ = 'lena';

% compressive imaging recontruction problem definition
subrate = 0.2; % subsampling rate
SNR_dB = [-5:5:30]; SNR_linear=10.^(SNR_dB./10);% input snr (the higher, the less noisy the problem is)  [-10:5:30]
opt.omega = 0; % snipe prior parameter (small values are generally better), see http://arxiv.org/abs/1312.3968 for details 

% scampi options
varNoiseDiffs = 1e-1; % initial noise variance associated to the dual variables
opt.learnNoise =1; % use the Bethe free energy to learn the noise variances
opt.dump_learn = 0.9; % damping for the learning of the noise variance (the larger, the higher the damping)
opt.dump_mes = 0.1; % damping of the AMP algorithm
opt.nb_iter = 300; % max number of iterations
opt.conv = 1e-5; % convergence criterion
opt.print = 10; % print infos every *opt.print* iterations
opt.showDynamics = 0; % show the reconstructed image in real time

% channel information
K = 1;
L = 3;
N1 = sqrt(N);
N2 = sqrt(N);
% D_y = 2.5*f;%equivalent lens length
% D_z = 2.5*f;%equivalent lens height
D_y = 12;%equivalent lens length
D_z = 12;%equivalent lens height
A = 20;% length aperture
V1 = 8;
V2 = 8;
%Num = 1;%实验次数
Q = 512;
%generate_matrix(N, N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Num=1;
for ii=1:length(SNR_dB)
  ii
    sigma2=1/SNR_linear(ii);
    error1 = 0; error2 = 0; error3 = 0;error4=0;
    
for iter = 1:Num;  
      iter
      %H_beam = beamspace_channel_lens(N1,N2,K,L,[],[],D_y,D_z,A);    % beamspace channel
      H_beam = test(N1,N2,K,L,[],[],D_y,D_z,A);    % beamspace channel
     for j = 1 : K
      h = reshape(H_beam,N,1) ;%channel vector for k th user
       original_image=reshape(h,N1,N2);
%       downscale = sqrt(N) / max(size(original_image) );
%       original_image = round(imresize(original_image, round(downscale * size(original_image) ) ) ); % downscale image if needed  
       [original_image1, mmin, mmax] = rescaleImage(original_image); % rescale image to have pixel values in [0,1]
      opt.signal = 255*original_image1(:);
   M = round(subrate * N); % number of lines of the original measurement matrix
if had
    % create the augmented operators for the cosparse TV analysis model with dual variables
    [op, opSq, opTr, opSqTr, l, ~, Ntot, rp, ~, rp2] = createAugmentedOp(N, subrate, 1);

    % compressive linear measurement    
    z = MultSeededHadamard(opt.signal, 1 / N, 1, 1, M, N, rp, [], rp2);          
else        
    % create the augmented operators for the cosparse TV analysis model with dual variables    
    F = randn(M, N) ./ sqrt(N); % random Gaussian operator
    [op, opSq, opTr, opSqTr, l, ~, Ntot] = createAugmentedOp(N, subrate, 1, F);

    % compressive linear measurement
    z = F * opt.signal;  
end  

% noise
%nvar = mean(z(z < max(abs(z) ) ).^2) * exp(-log(10) * SNR_dB(ii) / 10); % noise variance (avoiding to take into account the measurement associated to the first Hadamard mode)
%nvar = mean(z(z < max(abs(z) ) ).^2)/10^(SNR_dB(ii)/10);
%nvar
nvar=(norm(z)^2/length(z)/SNR_linear(ii));
noise = sqrt(nvar) * randn(M,1); % zero mean additive white Gaussian noise    
y = z + real(noise); 
yz = [y; sparse(l, 1) ]; % augmented measurement
 
%% scampi
% some fields required by scampi              
if opt.learnNoise; % noise variance (associated to pixels and dual variables)    
    opt.var_noise = [1e-1 * ones(M, 1); varNoiseDiffs * ones(l, 1) ];
else 
    opt.var_noise = [nvar * ones(M, 1); varNoiseDiffs * ones(l, 1) ];
end                
opt.M = numel(yz); % size of augmented measurement
opt.N = Ntot; % size of augmented signal
opt.Ntrue = N; % size of image
opt.Mtrue = M; % size of original measurement
[X_scampi, D_scampi] = scampi_solver(yz, op, opTr, opSq, opSqTr, opt,SNR_dB(ii)); % X_scampi is the reconstructed image, D_scampi, the estimation of the dual variables

%% ssd
Phi = (1/sqrt(Q))*(2*(rand(Q,N)<0.5) - 1);%             Phi = exp(1i*2*pi*rand(96,n));  % analog beamforming
%noise2 = sqrt(sigma2)*(randn(Q,1)+1i*randn(Q,1))/sqrt(2);
%Phi=F;
%nvar2 = mean(h(h < max(abs(h) ) ).^2) * exp(-log(10) * isnr(ii) / 10);
z1=Phi*h ;
nvar1=(norm(z1)^2/length(z1)/SNR_linear(ii));
noise1 = sqrt(nvar1) * randn(Q,1);
xx = z1+ noise1;

[h_hat2 n1 n2 support r] = SSD_simplify(xx,Phi,L,N1,N2,h,V1,V2);

%% ldamp
SamplingRate=subrate ;
AMP_iters=20;
VAMP_iters=20;
imsize=128; %512
n_DnCNN_layers=20;%Other option is 17
LoadNetworkWeights(n_DnCNN_layers);
T= round(SamplingRate * N);
F=(1/sqrt(T))*(2*(rand(T,N)<0.5) - 1);
z2=F*opt.signal;

B=@(x) F*x;
Bt=@(x) F.'*x;
height=N1;
width=N2;

MSEfxn = @(x_hat) MSE_cal(255*original_image1,reshape(x_hat,[height width]));
%nvar = mean(z(z < max(abs(z) ) ).^2) * exp(-log(10) * isnr(ii) / 10); % noise variance (avoiding to take into account the measurement associated to the first Hadamard mode)
nvar2=(norm(z2)^2/length(z2)/SNR_linear(ii));
noise2 = sqrt(nvar2) * randn(T,1); % zero mean additive white Gaussian noise    
y = z2 + real(noise2); 
[x_LDAMP,mmse_LDAMP(ii)] = DAMP_SNR(y,AMP_iters,N1,N2,'DnCNN',B,Bt,MSEfxn);
[x_LDAMP1,mmse_LDAMP1(ii)] = DAMP_SNR(y,AMP_iters,N1,N2,'BLS-GSM',B,Bt,MSEfxn); %Available options are NLM, Gauss, Bilateral, BLS-GSM, BM3D, fast-BM3D, BM3D-SAPCA, and DnCNN



error1 = error1 + norm(opt.signal - X_scampi,2)^2/norm(opt.signal,2)^2;
error2 = error2 + norm(h - h_hat2,2)^2/norm(h,2)^2;
error3=  error3 + norm(opt.signal - x_LDAMP(:),2)^2/norm(opt.signal,2)^2;
error4=  error4 + norm(opt.signal - x_LDAMP1(:),2)^2/norm(opt.signal,2)^2;

end
end
NMSE1(ii) = error1/K/Num;
NMSE2(ii) = error2/K/Num;
NMSE3(ii) = error3/K/Num;
NMSE4(ii) = error4/K/Num;

end
x_LDAMP_linear=reshape(x_LDAMP,N,1);
x_LDAMP1_linear=reshape(x_LDAMP1,N,1);
semilogy(SNR_dB,NMSE1,'b','Linewidth',1.5,'Displayname',['SCAMPI']);
hold on 
semilogy(SNR_dB,NMSE2,'r','Linewidth',1.5,'Displayname',['SD']);
hold on 
semilogy(SNR_dB,NMSE3,'k','Linewidth',1.5,'Displayname',['DnCNN']);
hold on 
semilogy(SNR_dB,NMSE4,'m','Linewidth',1.5,'Displayname',['NLM']);
hold on 

grid on
xlabel('SNR (dB)');
ylabel('NMSE (dB)');hold on

