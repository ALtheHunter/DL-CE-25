clear all 
clc
had = 1; % Do you want to use an Hadamard or a random measurement matrix?
f = 2;
N = (16*f)^2; % size of image (power of two with Hadamard operator)
% image_ = 'lena';

% compressive imaging recontruction problem definition
subrate = 0.1; % subsampling rate
isnr = 30;%[-20:10:30]; 
SNR_linear=10.^(isnr/10.);% input snr (the higher, the less noisy the problem is)  [-10:5:30]
opt.omega = 0; % snipe prior parameter (small values are generally better), see http://arxiv.org/abs/1312.3968 for details 

% scampi options
varNoiseDiffs = 1e-1; % initial noise variance associated to the dual variables
opt.learnNoise = 1; % use the Bethe free energy to learn the noise variances
opt.dump_learn = 0.9; % damping for the learning of the noise variance (the larger, the higher the damping)
opt.dump_mes = 0.1; % damping of the AMP algorithm
opt.nb_iter = 300; % max number of iterations
opt.conv = 1e-5; % convergence criterion
opt.print = 10; % print infos every *opt.print* iterations
opt.showDynamics = 0; % show the reconstructed image in real time

% channel information
K = 4;
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
N_iter = 1;
Q = 512;
%generate_matrix(N, N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i_snr=1:length(isnr)
     i_snr
    sigma2=1/SNR_linear(i_snr);
    error1 = 0; error2 = 0; error3 = 0;
   
    
    for iter = 1:N_iter;  
      iter
      %H_beam = beamspace_channel_lens(N1,N2,K,L,[],[],D_y,D_z,A);    % beamspace channel
      H_beam = test(N1,N2,K,L,[],[],D_y,D_z,A);    % beamspace channel
     for j = 1 : K
      h = H_beam(:,j);%channel vector for k th user
       original_image=reshape(h,N1,N2);
%       downscale = sqrt(N) / max(size(original_image) );
%       original_image = round(imresize(original_image, round(downscale * size(original_image) ) ) ); % downscale image if needed  
      [original_image1, mmin, mmax] = rescaleImage(original_image); % rescale image to have pixel values in [0,1]
      opt.signal = original_image1(:);


      
M = round(subrate * N); % number of lines of the original measurement matrix
if had
    % create the augmented operators for the cosparse TV analysis model with dual variables
    [op, opSq, opTr, opSqTr, l, ~, Ntot, rp, ~, rp2] = createAugmentedOp(N, subrate, 1);
    index = rp;
    W_r = recover_hadamard (rp, N, M); %recover hadamard matrix from fast hadamard algorithm
    % compressive linear measurement    
    yt = MultSeededHadamard(opt.signal, 1 / N, 1, 1, M, N, rp, [], rp2);
    input_signal(rp2) = opt.signal;
    yr = W_r * input_signal'; 
else        
    % create the augmented operators for the cosparse TV analysis model with dual variables    
    F = randn(M, N) ./ sqrt(N); % random Gaussian operator
    [op, opSq, opTr, opSqTr, l, ~, Ntot] = createAugmentedOp(N, subrate, 1, F);

    % compressive linear measurement
    y = F * opt.signal;  
end  

% noise
nvar = mean(yr(yr < max(abs(yr) ) ).^2) * exp(-log(10) * isnr(i_snr) / 10); % noise variance (avoiding to take into account the measurement associated to the first Hadamard mode)
nvar
noise = sqrt(nvar) * randn(M,1); % zero mean additive white Gaussian noise    
y = yr + real(noise); 
%yz = [y; sparse(l, 1) ]; 
yz = [y; zeros(l,1)];  % augmented measurement

% augmented measurement matrix 
W_aug = augment_matrix(N1,N2,N,M,l,W_r);

% ssd
Phi = (1/sqrt(Q))*(2*(rand(Q,N)<0.5) - 1);%             Phi = exp(1i*2*pi*rand(96,n));  % analog beamforming
%noise2 = sqrt(sigma2)*(randn(Q,1)+1i*randn(Q,1))/sqrt(2);
nvar2 = mean(h(h < max(abs(h) ) ).^2) * exp(-log(10) * isnr(i_snr) / 10);
noise2 = sqrt(nvar2) * randn(Q,1);
xx = Phi*h + noise2;
[h_hat2 n1 n2 support r] = SSD_simplify(xx,Phi,L,N1,N2,h,V1,V2);

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

% scoampi 
sigPow = (mean(opt.signal))^2;
[X_scampi, D_scampi,prior] = oamp_solver(yz, op, opTr, opSq, opSqTr, opt,isnr(i_snr),W_aug,sigPow,N,l,nvar); % X_scampi is the reconstructed image, D_scampi, the estimation of the dual variables
nsnr1(i_snr) = 10 * log10(norm(opt.signal(:) - X_scampi(1 : N) ).^2 ./ norm(opt.signal).^2); % final nsnr
disp(sprintf('final nsnr1=%0.2e', nsnr1) );
% subplot(1, 2, 1); imshow(original_image); title('original', 'Fontsize', 15);hold on
% subplot(1, 2, 2); imshow(reshape(X_scampi, sqrt(N), sqrt(N) ) ); title('scampi estimate', 'Fontsize', 15);hold on
error1 = error1 + norm(opt.signal - X_scampi,2)^2/norm(opt.signal,2)^2;
error2 = error2 + norm(h - h_hat2,2)^2/norm(h,2)^2;
nsnr2(i_snr) = 10 * log10(norm(h(:) - h_hat2(1 : N) ).^2 ./ norm(h).^2); % final nsnr
%recoveredimage = recoverimage(X_scampi, mmin, mmax);
%hest = reshape(recoveredimage,N1,N2);
%error3 = error3 + norm(original_image(:) - recoveredimage,2)^2/norm(original_image(:),2)^2;
     end
   end
   
NMSE1(i_snr) = error1/K/N_iter;
NMSE2(i_snr) = error2/K/N_iter;
%NMSE3(i_snr) = error3/K/N_iter;
end

figure(1)
semilogy(isnr,NMSE1,'b','Linewidth',1.5);
hold on 
semilogy(isnr,NMSE2,'r','Linewidth',1.5);
%hold on 
%semilogy(isnr,NMSE3,'y','Linewidth',1.5);
grid on
xlabel('SNR (dB)');
ylabel('NMSE (dB)');

% figure(3)
% plot(isnr,nsnr1,'b','Linewidth',1.5);
% hold on 
% plot(isnr,nsnr2,'r','Linewidth',1.5);
% grid on
% xlabel('ISNR (dB)');
% ylabel('NSNR (dB)');