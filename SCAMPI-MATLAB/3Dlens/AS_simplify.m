function [h_hat,support] = AS_simplify(x,Phi,n1,n2,N1,N2,V1,V2)
n = N1*N2;
N=N1;
for i = -V1/2:1:V1/2
    i1 = i + V1/2 + 1;
    for j = -V2/2:1:V2/2
        j1 = j + V2/2 +1;
        select((i1-1)*(V2+1) + j1,1) = n1 + i;
        select((i1-1)*(V2+1) + j1,2) = n2 + j;
    end
end
for i = 1:length(select)
    if select(i,1) > N1
       select(i,1) = select(i,1) - N1;
    elseif select(i,1)<1
       select(i,1) = select(i,1) + N1;
    end
end
for i = 1:length(select)
    if select(i,2) > N2
       select(i,2) = select(i,2) - N2;
    elseif select(i,2)<1
       select(i,2) = select(i,2) + N2;
    end
end
for i = 1:length(select)
    support(i) = N*(select(i,2)-1) + select(i,1);
end
Phi2 = Phi(:,unique(support));
h_hat2 = inv(Phi2'*Phi2)*Phi2'*x;
temp = zeros(n,1);
temp(unique(support)) = h_hat2;
b1 = n1 - V1/2;
if b1 > N1
   b1 = b1 - N1;
elseif b1 < 1
   b1 = b1 + N1;
end
b2 = n1 + V1/2;
if b2 > N1
   b2 = b2 - N1;
elseif b2 < 1
   b2 = b2 + N1;
end
b3 = n2 - V2/2;
if b3 > N2
   b3 = b3 - N2;
elseif b3 < 1
   b3 = b3 + N2;
end 
b4 = n2 + V2/2;
if b4 > N2
   b4 = b4 - N2;
elseif b4 < 1
   b4 = b4 + N2;
end
h_hat = temp;