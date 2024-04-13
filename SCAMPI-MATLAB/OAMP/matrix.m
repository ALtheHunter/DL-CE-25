function y = matrix(q,n,m)
y = zeros(q,n);
for i = 1:q
   y(i,i)=1;
   y(i,m+i)=-1;
end
end