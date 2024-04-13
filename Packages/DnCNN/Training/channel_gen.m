%%% This script generates the channel data for 3D lens antenna array-based mmWave system
%%% The aim is to produce training data for DnCNN
clear all
clc
NumTraining_data=16640;% the number of data
Numvali_data=6400;% the number of data
K = 1;
L = 4;
M = 64;% The number of antennas
N = 64;% The number of antennas
Training_h=zeros(M,N,1,NumTraining_data);
h_norm=zeros(M,N,1,NumTraining_data);
vali_h=zeros(M,N,1,Numvali_data);
h_norm_vali=zeros(M,N,1,Numvali_data);
% D_y = 2.5*f;%equivalent lens length
% D_z = 2.5*f;%equivalent lens height
D_y = 12;%equivalent lens length
D_z = 12;%equivalent lens height
A = 20;% length aperture
%N_iter = 5;
%Q = 512;
% SNR=50dB;
% noise= 1/10^(SNR/10);


%generate_matrix(N, N);
for ii=1:NumTraining_data
  Training_h(:,:,1,ii)= test(M,N,K,L,[],[],D_y,D_z,A) ;
  [h_norm(:,:,1,ii)] = rescaleImage(Training_h(:,:,1,ii));
end
set=uint8(ones(1,NumTraining_data));

inputs=single(h_norm);


save('D:\LDAMP_for_Rice\D-AMP_Toolbox-master\Packages\DnCNN\Training\Training_h.mat','inputs','set');


for ii=1: Numvali_data
    
  vali_h(:,:,1,ii)= test(M,N,K,L,[],[],D_y,D_z,A);
  [h_norm_vali(:,:,1,ii)] = rescaleImage(vali_h(:,:,1,ii));
  
end


inputs=single(h_norm_vali);

% vali_h=[h_norm(:,:,1,1:Numvali_data)];
 set=uint8(zeros(1,Numvali_data));

save('D:\LDAMP_for_Rice\D-AMP_Toolbox-master\Packages\DnCNN\Training\vali_h','inputs','set');