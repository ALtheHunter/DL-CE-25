%%% This script generates the channel data for 3D lens antenna array-based mmWave system
%%% The aim is to produce training data for DnCNN
clear all
clc
NumTraining_data=16640;% the number of data
Numvali_data=6528;% the number of data
K = 1;
L = 3;
M = 50;% The number of antennas
N = 50;% The number of antennas
Training_h=zeros(M,N,1,NumTraining_data);
h_norm=zeros(M,N,1,NumTraining_data);
vali_h=zeros(M,N,1,Numvali_data);
% D_y = 2.5*f;%equivalent lens length
% D_z = 2.5*f;%equivalent lens height
D_y = 12;%equivalent lens length
D_z = 12;%equivalent lens height
A = 20;% length aperture
%N_iter = 5;
%Q = 512;
% SNR=50dB;
% noise= 1/10^(SNR/10);
set=ones(1,NumTraining_data);

%generate_matrix(N, N);
for ii=1:NumTraining_data
  Training_h(:,:,1,ii)= test(M,N,K,L,[],[],D_y,D_z,A);
  [h_norm(:,:,1,ii)] = rescaleImage(Training_h(:,:,1,ii)); % rescale image to have pixel values in [0,1]
end


inputs=single(h_norm);


save('C:\Users\dell\Desktop\OAMP_for_channel_estimation\SCAMPI-MATLAB\scampi-vs-ssd\Training_h','inputs','set');


for ii=1: Numvali_data
  vali_h(:,:,1,ii)= test(M,N,K,L,[],[],D_y,D_z,A) ;
end


inputs=single(vali_h);
set=zeros(1,Numvali_data);



save('C:\Users\dell\Desktop\OAMP_for_channel_estimation\SCAMPI-MATLAB\scampi-vs-ssd\vali_h','inputs','set');