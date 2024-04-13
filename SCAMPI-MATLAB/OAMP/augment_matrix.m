function [W_aug, D] = augment_matrix(N1,N2,N,M,l,W_r)
D_v = matrix(N1*(N2-1), N, 1); 
D_hh = matrix(N2*(N1-1), N, M);
D_hh = D_hh(1:N-M,1:N);
D_h = zeros(N1*(N2-1),N);
D_h(1:N-M,1:N) = D_hh;
D = zeros(l,N);
D(1:N1*(N2-1),1:N) = D_h;
D(N1*(N2-1)+1:N1*(N2-1)+N2*(N1-1),1:N) = D_v;
I = -eye(l);
O = zeros(M,l);
W_r_aug = [W_r,O];
D_aug = [D,I];
% for j=1:N+l; No(j)=norm(D_aug(:,j)); end 
% D_aug = D_aug*diag(1./No);    
W_aug = zeros(M+l,N+l);
W_aug(1:M,1:N+l) = W_r_aug;
W_aug(1+M:M+l,1:N+l) = D_aug;