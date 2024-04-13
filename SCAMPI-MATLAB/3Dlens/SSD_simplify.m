function [h_out n1 n2 support r] = SSD_simplify(x,Phi,L,N1,N2,h,V1,V2)
n = N1*N2;
x_temp = x;
support=[];
for l = 1 : L
    r = Phi'*x_temp;
    z = reshape(r,N1,N2);
    [n1 n2] = find(z==max(max(z)));

    %%% adaptive support detection
    [temp,support(l,:)] = AS_simplify(x,Phi,n1,n2,N1,N2,V1,V2);
 
%     for i = -V2/2:1:V2/2-1
%         i1 = i + V2/2 + 1;
%         for j = -V1/2:1:V1/2-1
%             j1 = j + V1/2 +1;
%             select((i1-1)*V1 + j1,1) = n1 + i;
%             select((i1-1)*V1 + j1,2) = n2 + j;
%         end
%     end
%     
%     
% %     for i = -V1/2:1:V1/2
% %         i1 = i + V1/2 + 1;
% %         for j = -V2/2:1:V2/2
% %             j1 = j + V2/2 +1;
% %             select((i1-1)*(V1+1) + j1,1) = n1 + i;
% %             select((i1-1)*(V1+1) + j1,2) = n2 + j;
% %         end
% %     end
% %     select =
% %     [n1-2,n2;n1-1,n2-1;n1-1,n2;n1-1,n2+1;n1,n2-2;n1,n2-1;n1,n2;n1,n2+1;n1,n2+2;n1+1,n2-1;n1+1,n2;n1+1,n2+1;n1+2,n2];
%     for i = 1:length(select)
%         if select(i,1) > N1
%            select(i,1) = select(i,1) - N1;
%         elseif select(i,1)<1
%            select(i,1) = select(i,1) + N1;
%         end
%     end
%     for i = 1:length(select)
%         if select(i,2) > N2
%            select(i,2) = select(i,2) - N2;
%         elseif select(i,2)<1
%            select(i,2) = select(i,2) + N2;
%         end
%     end
%     for i = 1:length(select)
%         support(l,i) = 32*(select(i,2)-1) + select(i,1);
%     end
%     Phi2 = Phi(:,unique(support(l,:)));
%     h_hat2 = inv(Phi2'*Phi2)*Phi2'*x_temp;
%     temp = zeros(n,1);
%     temp(unique(support(l,:))) = h_hat2;
    if l>=1
       x_temp = x_temp - Phi*temp;
    end
end
%h_out = temp;

[g o]=size(support);
 support_tot = unique(reshape(support,g*o,1));
 Phi_final = Phi(:,support_tot);
 est =  inv(Phi_final'*Phi_final)*Phi_final'*x;
% % est = inv(Phi_final'*Phi_final+sigma2*eye(length(select_final)))*Phi_final'*x;
 h_out = zeros(n,1);
 h_out(support_tot) = est;

