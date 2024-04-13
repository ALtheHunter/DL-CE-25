clear all 
clc
had = 1; % Do you want to use an Hadamard or a random measurement matrix?
f = 8;
N = (16*f)^2; % size of image (power of two with Hadamard operator)
% image_ = 'lena';

% compressive imaging recontruction problem definition
subrate = 0.1; % subsampling rate
isnr = [-20:10:30]; SNR_linear=10.^(isnr/10.);% input snr (the higher, the less noisy the problem is)  [-10:5:30]
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
L = 1;
N1 = sqrt(N);
N2 = sqrt(N);
D_y = 2.5*f;%equivalent lens length
D_z = 2.5*f;%equivalent lens height
A = 20;% length aperture
V1 = 8;
V2 = 8;
N_iter=5;
Q=8192;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i_snr=1:length(isnr)
     i_snr
    sigma2=1/SNR_linear(i_snr);
    error1 = 0; error2 = 0; 
   for iter = 1:N_iter;  
      iter
      H_beam = test(N1,N2,K,L,[],[],D_y,D_z,A);    % beamspace channel
     for j = 1 : K
      h = H_beam(:,j);%channel vector for k th user
% ssd
Phi = (1/sqrt(Q))*(2*(rand(Q,N)<0.5) - 1);%           
nvar2 = mean(h(h < max(abs(h) ) ).^2) * exp(-log(10) * isnr(i_snr) / 10);
noise2 = sqrt(nvar2) * randn(Q,1);
xx = Phi*h + noise2;
[h_hat2 n1 n2 support r] = SSD_simplify(xx,Phi,L,N1,N2,h,V1,V2);
error2 = error2 + norm(h - h_hat2,2)^2/norm(h,2)^2;
nsnr2(i_snr) = 10 * log10(norm(h(:) - h_hat2(1 : N) ).^2 ./ norm(h).^2); % final nsnr

     end
   end
   

NMSE2(i_snr) = error2/K/N_iter;
end

figure(1)
semilogy(isnr,NMSE2,'r','Linewidth',1.5);
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